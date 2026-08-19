#!/usr/bin/env python3
from pathlib import Path
import re
import sys

root = Path(sys.argv[1]).resolve()

# Build against the official NeoForge 26.2 MDK baseline.  Keep the first pass
# deliberately dependency-light: optional integrations are disabled until the
# core is source-compatible with 26.2.
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
        content {
            includeGroup("com.github.FiguraMC.luaj")
        }
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
        create("tacz") {
            sourceSet(sourceSets["main"])
        }
    }
}

sourceSets {
    main {
        java {
            // Re-enable integrations one at a time after the core compiles.
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
}
'''
(root / "build.gradle.kts").write_text(build, encoding="utf-8")

# 26.2 official MDK currently uses Gradle 9.2.1 and Java 25.
wrapper = root / "gradle/wrapper/gradle-wrapper.properties"
wrapper.write_text("""distributionBase=GRADLE_USER_HOME\ndistributionPath=wrapper/dists\ndistributionUrl=https\\://services.gradle.org/distributions/gradle-9.2.1-bin.zip\nnetworkTimeout=10000\nvalidateDistributionUrl=true\nzipStoreBase=GRADLE_USER_HOME\nzipStorePath=wrapper/dists\n""", encoding="utf-8")

# Remove old 1.21.1 build metadata that can confuse Gradle during migration.
for p in [root / "gradle/libs.versions.toml"]:
    if p.exists():
        p.rename(p.with_suffix(p.suffix + ".disabled"))

# The optional compatibility registry directly references integrations that are
# intentionally excluded during the core migration pass. Replace it with a
# stable no-op bridge while preserving its public init hooks where possible.
compat = root / "src/main/java/com/tacz/guns/init/CompatRegistry.java"
if compat.exists():
    text = compat.read_text(encoding="utf-8")
    # Discover the package/class and all public static void zero-arg methods,
    # then preserve those signatures as no-ops to avoid touching callers.
    methods = re.findall(r'public\s+static\s+void\s+(\w+)\s*\(\s*\)', text)
    methods = list(dict.fromkeys(methods))
    body = ["package com.tacz.guns.init;", "", "/**", " * 26.2 port bridge. Optional mod integrations are registered after the core", " * migration is complete; keeping these hooks no-op preserves call sites.", " */", "public final class CompatRegistry {", "    private CompatRegistry() {}"]
    if not methods:
        methods = ["init"]
    for m in methods:
        body += [f"    public static void {m}() {{}}"]
    body += ["}", ""]
    compat.write_text("\n".join(body), encoding="utf-8")

# Disable KubeJS service discovery while the package is excluded.
for rel in ["src/main/resources/kubejs.plugins.txt"]:
    p = root / rel
    if p.exists():
        p.write_text("", encoding="utf-8")

print("Applied TACZ NeoForge 26.2 migration pass 1")
