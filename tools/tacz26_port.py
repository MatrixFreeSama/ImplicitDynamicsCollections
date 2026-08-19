#!/usr/bin/env python3
from pathlib import Path
import re
import shutil
import sys

root = Path(sys.argv[1]).resolve()
fabric = root.parent / "fabric26"

# Build against the official NeoForge 26.2 MDK baseline. Optional third-party
# integrations stay disabled until the core mod is source-compatible.
build = r'''plugins {
    java
    id("net.neoforged.moddev") version "2.0.144"
}

group = "com.tacz"
version = "1.1.8-neoforge-26.2-port"

base {
    archivesName.set("tacz-neoforge-26.2")
}

java {
    toolchain.languageVersion = JavaLanguageVersion.of(25)
    withSourcesJar()
}

repositories {
    mavenCentral()
    maven("https://jitpack.io") {
        content { includeGroup("com.github.FiguraMC.luaj") }
    }
}

neoForge {
    version = "26.2.0.59"

    runs {
        create("client") {
            client()
            gameDirectory = file("run/client")
        }
        create("server") {
            server()
            gameDirectory = file("run/server")
            programArgument("--nogui")
        }
    }

    mods {
        create("tacz") { sourceSet(sourceSets["main"]) }
    }
}

sourceSets {
    main {
        java {
            // Integrations are optional and are re-enabled after the core port.
            exclude("com/tacz/guns/compat/**")
        }
    }
}

dependencies {
    implementation("org.apache.commons:commons-math3:3.6.1")
    implementation("com.github.FiguraMC.luaj:luaj-core:3.0.8-figura")
    implementation("com.github.FiguraMC.luaj:luaj-jse:3.0.8-figura")
    implementation("org.apache.bcel:bcel:6.6.1")
}

tasks.withType<JavaCompile>().configureEach {
    options.encoding = "UTF-8"
    options.release.set(25)
    options.compilerArgs.addAll(listOf("-Xmaxerrs", "2000", "-Xmaxwarns", "2000"))
}
'''
(root / "build.gradle.kts").write_text(build, encoding="utf-8")

wrapper = root / "gradle/wrapper/gradle-wrapper.properties"
wrapper.write_text("""distributionBase=GRADLE_USER_HOME\ndistributionPath=wrapper/dists\ndistributionUrl=https\\://services.gradle.org/distributions/gradle-9.2.1-bin.zip\nnetworkTimeout=10000\nvalidateDistributionUrl=true\nzipStoreBase=GRADLE_USER_HOME\nzipStorePath=wrapper/dists\n""", encoding="utf-8")

libs = root / "gradle/libs.versions.toml"
if libs.exists():
    libs.rename(libs.with_suffix(libs.suffix + ".disabled"))

# Pull verified vanilla/API migrations from the working 26.2 Fabric port, but
# only for files that are platform-neutral on both sides. This avoids guessing
# hundreds of Mojang 26.x rendering/data API changes while keeping NeoForge's
# registrations, networking and event wiring intact.
if fabric.exists():
    neo_java = root / "src/main/java"
    fab_java = fabric / "src/main/java"
    neo_files = {p.relative_to(neo_java).as_posix(): p for p in neo_java.rglob("*.java")}
    fab_files = {p.relative_to(fab_java).as_posix(): p for p in fab_java.rglob("*.java")}
    copied = 0
    for rel in sorted(set(neo_files) & set(fab_files)):
        old = neo_files[rel].read_text(encoding="utf-8", errors="ignore")
        new = fab_files[rel].read_text(encoding="utf-8", errors="ignore")
        fabric_imports = []
        for line in new.splitlines():
            s = line.strip()
            if s.startswith("import net.fabricmc") or s.startswith("import cn.sh1rocu"):
                fabric_imports.append(s.removeprefix("import ").removesuffix(";"))
        only_env_annotations = all(x in {
            "net.fabricmc.api.EnvType", "net.fabricmc.api.Environment"
        } for x in fabric_imports)
        platform_area = any(seg in rel for seg in (
            "/init/", "/network/", "/event/", "/mixin/", "/compat/"
        ))
        if only_env_annotations and "net.neoforged" not in old and not platform_area:
            new = new.replace("import net.fabricmc.api.EnvType;\n", "")
            new = new.replace("import net.fabricmc.api.Environment;\n", "")
            new = re.sub(r"\s*@Environment\(EnvType\.(?:CLIENT|SERVER)\)\s*", "\n", new)
            neo_files[rel].write_text(new, encoding="utf-8")
            copied += 1
    print(f"Copied {copied} platform-neutral 26.2 source files")

# Mechanical Mojang 26.2 renames for the NeoForge-specific files retained from
# the 1.21.1 port. Identifier replaced ResourceLocation in 26.2; Util moved to
# net.minecraft.util; RenderType moved under renderer.rendertype.
for p in (root / "src/main/java").rglob("*.java"):
    text = p.read_text(encoding="utf-8", errors="ignore")
    text = text.replace("net.minecraft.resources.ResourceLocation", "net.minecraft.resources.Identifier")
    text = re.sub(r"\bResourceLocation\b", "Identifier", text)
    text = text.replace("import net.minecraft.Util;", "import net.minecraft.util.Util;")
    text = text.replace("import net.minecraft.client.renderer.RenderType;", "import net.minecraft.client.renderer.rendertype.RenderType;")
    p.write_text(text, encoding="utf-8")

# KubeJS is optional. Keep the public event-poster surface but make it a no-op
# while the KubeJS compatibility package is intentionally excluded.
kube = root / "src/main/java/com/tacz/guns/api/event/common/KubeJSGunEventPoster.java"
if kube.exists():
    kube.write_text('''package com.tacz.guns.api.event.common;\n\nimport net.neoforged.bus.api.Event;\n\npublic interface KubeJSGunEventPoster<E extends Event> {\n    default void postEventToKubeJS(E event) {}\n    default void postClientEventToKubeJS(E event) {}\n    default void postServerEventToKubeJS(E event) {}\n}\n''', encoding="utf-8")

# Optional integration registry bridge. Preserve the zero-argument hooks used by
# core setup without linking excluded compatibility implementations.
compat = root / "src/main/java/com/tacz/guns/init/CompatRegistry.java"
if compat.exists():
    text = compat.read_text(encoding="utf-8", errors="ignore")
    methods = list(dict.fromkeys(re.findall(r'public\s+static\s+void\s+(\w+)\s*\(\s*\)', text)))
    if not methods:
        methods = ["init"]
    body = [
        "package com.tacz.guns.init;", "",
        "public final class CompatRegistry {", "    private CompatRegistry() {}"
    ]
    for m in methods:
        body.append(f"    public static void {m}() {{}}")
    body += ["}", ""]
    compat.write_text("\n".join(body), encoding="utf-8")

for rel in ["src/main/resources/kubejs.plugins.txt"]:
    p = root / rel
    if p.exists():
        p.write_text("", encoding="utf-8")

print("Applied TACZ NeoForge 26.2 migration pass 2")
