#!/usr/bin/env python3
from pathlib import Path
import re
import shutil
import sys

root = Path(sys.argv[1]).resolve()
fabric = root.parent / "fabric26"
java_root = root / "src/main/java"

build = r'''plugins {
    java
    id("net.neoforged.moddev") version "2.0.144"
}

group = "com.tacz"
version = "1.1.8-neoforge-26.2-port"

base { archivesName.set("tacz-neoforge-26.2") }

java {
    toolchain.languageVersion = JavaLanguageVersion.of(25)
    withSourcesJar()
}

repositories {
    mavenCentral()
    maven("https://jitpack.io") { content { includeGroup("com.github.FiguraMC.luaj") } }
}

neoForge {
    version = "26.2.0.59"
    runs {
        create("client") { client(); gameDirectory = file("run/client") }
        create("server") { server(); gameDirectory = file("run/server"); programArgument("--nogui") }
    }
    mods { create("tacz") { sourceSet(sourceSets["main"]) } }
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
    options.compilerArgs.addAll(listOf("-Xmaxerrs", "3000", "-Xmaxwarns", "3000"))
}
'''
(root / "build.gradle.kts").write_text(build, encoding="utf-8")

(root / "gradle/wrapper/gradle-wrapper.properties").write_text(
    "distributionBase=GRADLE_USER_HOME\n"
    "distributionPath=wrapper/dists\n"
    "distributionUrl=https\\://services.gradle.org/distributions/gradle-9.2.1-bin.zip\n"
    "networkTimeout=10000\nvalidateDistributionUrl=true\n"
    "zipStoreBase=GRADLE_USER_HOME\nzipStorePath=wrapper/dists\n", encoding="utf-8")
libs = root / "gradle/libs.versions.toml"
if libs.exists(): libs.rename(libs.with_suffix(libs.suffix + ".disabled"))

def strip_env(text: str) -> str:
    text = text.replace("import net.fabricmc.api.EnvType;\n", "")
    text = text.replace("import net.fabricmc.api.Environment;\n", "")
    text = re.sub(r"\s*@Environment\(EnvType\.(?:CLIENT|SERVER)\)\s*", "\n", text)
    return text

def imports_with_prefix(text: str, prefixes):
    out=[]
    for line in text.splitlines():
        s=line.strip()
        if s.startswith("import "):
            imp=s.removeprefix("import ").removesuffix(";")
            if imp.startswith(prefixes): out.append(imp)
    return out

# Use the already-working 26.2 Fabric port as the source of truth for Mojang's
# 26.2 API changes, while retaining NeoForge-specific event/registry/network files.
if fabric.exists():
    fab_root = fabric / "src/main/java"
    neo = {p.relative_to(java_root).as_posix(): p for p in java_root.rglob("*.java")}
    fab = {p.relative_to(fab_root).as_posix(): p for p in fab_root.rglob("*.java")}
    copied_common=0
    for rel in sorted(set(neo) & set(fab)):
        old=neo[rel].read_text(encoding="utf-8", errors="ignore")
        new=fab[rel].read_text(encoding="utf-8", errors="ignore")
        fi=imports_with_prefix(new,("net.fabricmc","cn.sh1rocu"))
        ni=imports_with_prefix(old,("net.neoforged",))
        fabric_is_neutral=all(x in {"net.fabricmc.api.EnvType","net.fabricmc.api.Environment"} for x in fi)
        neo_is_neutral=all(x in {"net.neoforged.api.distmarker.Dist","net.neoforged.api.distmarker.OnlyIn"} for x in ni)
        if fabric_is_neutral and neo_is_neutral:
            neo[rel].write_text(strip_env(new), encoding="utf-8")
            copied_common += 1

    # New 26.2 helper/render/data classes which do not exist in the 1.21.1
    # NeoForge tree. Only loader-neutral com.tacz classes are admitted.
    copied_new=0
    allowed_prefixes=("java.","javax.","net.minecraft.","net.fabricmc.api.","com.tacz.",
        "org.jetbrains.","org.joml.","com.mojang.","cn.sh1rocu.","com.google.",
        "org.apache.","org.slf4j.","org.luaj.","org.lwjgl.","io.netty.","it.unimi.dsi.","static ")
    for rel in sorted(set(fab)-set(neo)):
        if not rel.startswith("com/tacz/guns/") or "/compat/" in rel or "/mixin/" in rel:
            continue
        text=fab[rel].read_text(encoding="utf-8", errors="ignore")
        fi=imports_with_prefix(text,("net.fabricmc","cn.sh1rocu"))
        if not all(x in {"net.fabricmc.api.EnvType","net.fabricmc.api.Environment"} for x in fi):
            continue
        bad=False
        for line in text.splitlines():
            s=line.strip()
            if not s.startswith("import "): continue
            imp=s.removeprefix("import ").removesuffix(";")
            if not imp.startswith(allowed_prefixes):
                bad=True; break
        if bad: continue
        dst=java_root/rel
        dst.parent.mkdir(parents=True,exist_ok=True)
        dst.write_text(strip_env(text),encoding="utf-8")
        copied_new += 1
    print(f"Copied {copied_common} common + {copied_new} new loader-neutral 26.2 Java files")

# Optional integrations are not allowed to hold the core port hostage. Remove
# their old 1.21.1 implementations and optional mixins, then provide small 26.2
# facades. Iris remains a real reflection-based compatibility bridge.
compat_dir=java_root/"com/tacz/guns/compat"
if compat_dir.exists(): shutil.rmtree(compat_dir)
for rel in [
    "com/tacz/guns/mixin/client/ar", "com/tacz/guns/mixin/client/iris",
    "com/tacz/guns/mixin/carryon", "com/tacz/guns/mixin/compat"
]:
    p=java_root/rel
    if p.exists(): shutil.rmtree(p)

def put(rel, text):
    p=java_root/rel; p.parent.mkdir(parents=True,exist_ok=True); p.write_text(text,encoding="utf-8")

put("com/tacz/guns/compat/ar/ARCompat.java", '''package com.tacz.guns.compat.ar;\nimport com.mojang.blaze3d.vertex.PoseStack;\nimport com.mojang.blaze3d.vertex.VertexConsumer;\npublic final class ARCompat { public static boolean LOADED=false; public static void init(){} public static boolean shouldAccelerate(){return false;} public static boolean isAccelerated(VertexConsumer v){return false;} public static void setRenderingLevel(){} public static void resetRenderingLevel(){} public static void setRenderLayer(int l){} public static void setRenderBeforeFunction(Runnable r){} public static void setRenderAfterFunction(Runnable r){} public static void resetRenderLayer(){} public static void resetRenderBeforeFunction(){} public static void resetRenderAfterFunction(){} public static void disableAcceleration(){} public static void resetAcceleration(){} public static void renderLaser(VertexConsumer v,float z,float w,boolean f,PoseStack p,int c){} }\n''')
put("com/tacz/guns/compat/controllable/ControllableCompat.java", '''package com.tacz.guns.compat.controllable;\nimport com.tacz.guns.api.item.gun.FireMode; import net.minecraft.world.item.ItemStack;\npublic final class ControllableCompat { public static void init(){} public static void onGunShoot(ItemStack s, FireMode m){} }\n''')
put("com/tacz/guns/compat/cloth/MenuIntegration.java", '''package com.tacz.guns.compat.cloth;\nimport net.minecraft.client.gui.screens.Screen;\npublic final class MenuIntegration { public static Screen getConfigScreen(Screen parent){ return null; } }\n''')
put("com/tacz/guns/compat/firstperson/FirstPersonAnimationCompat.java", '''package com.tacz.guns.compat.firstperson;\nimport net.minecraft.client.player.LocalPlayer; import net.minecraft.world.item.ItemStack;\npublic final class FirstPersonAnimationCompat { public static void init(){} public static ItemStack getMainRenderStack(LocalPlayer p){return p.getMainHandItem();} public static boolean isTaczViewmodel(ItemStack s){return true;} public static void beginDirectArmRender(){} public static void endDirectArmRender(){} }\n''')
put("com/tacz/guns/compat/immediatelyfast/ImmediatelyFastCompat.java", '''package com.tacz.guns.compat.immediatelyfast;\nimport net.minecraft.world.item.ItemStack; public final class ImmediatelyFastCompat { public static void init(){} public static void renderHotbarItem(ItemStack s,boolean pre){} public static boolean isInstalled(){return false;} }\n''')
put("com/tacz/guns/compat/shouldersurfing/ShoulderSurfingCompat.java", '''package com.tacz.guns.compat.shouldersurfing; public final class ShoulderSurfingCompat { public static void init(){} public static boolean showCrosshair(){return false;} public static boolean isInstalled(){return false;} }\n''')
put("com/tacz/guns/compat/zoomify/ZoomifyCompat.java", '''package com.tacz.guns.compat.zoomify; public final class ZoomifyCompat { public static void init(){} public static double getFov(double f,float t){return f;} }\n''')
put("com/tacz/guns/compat/playeranimator/PlayerAnimatorCompat.java", '''package com.tacz.guns.compat.playeranimator;\nimport com.tacz.guns.client.resource.GunDisplayInstance; import net.minecraft.world.entity.LivingEntity;\npublic final class PlayerAnimatorCompat { public static void init(){} public static boolean isInstalled(){return false;} public static boolean hasPlayerAnimator3rd(LivingEntity e,GunDisplayInstance d){return false;} public static void stopAllAnimation(LivingEntity e){} public static void stopAllAnimation(LivingEntity e,int f){} public static void playAnimation(LivingEntity e,GunDisplayInstance d,float l){} public static void registerReloadListener(Object o){} }\n''')
put("com/tacz/guns/compat/carryon/CarryOnReflection.java", '''package com.tacz.guns.compat.carryon;\nimport net.minecraft.core.BlockPos; import net.minecraft.core.HolderLookup; import net.minecraft.world.entity.player.Player; import net.minecraft.world.level.block.entity.BlockEntity; import net.minecraft.world.level.block.state.BlockState;\npublic final class CarryOnReflection { public static BlockState getCarriedBlock(Player p){return null;} public static BlockEntity getCarriedBlockEntity(Player p,BlockPos b,HolderLookup.Provider h){return null;} }\n''')

put("com/tacz/guns/compat/iris/IrisCompat.java", '''package com.tacz.guns.compat.iris;\nimport com.mojang.blaze3d.pipeline.RenderPipeline; import net.minecraft.client.renderer.SubmitNodeCollector; import net.neoforged.fml.ModList;\npublic final class IrisCompat {\n private static boolean iris(){ try{return ModList.get().isLoaded("iris");}catch(Throwable t){return false;} }\n public static void initCompat(){} public static boolean isRenderShadow(){ if(!iris())return false; try{Class<?> c=Class.forName("net.irisshaders.iris.shadows.ShadowRenderingState"); return (Boolean)c.getMethod("areShadowsCurrentlyBeingRendered").invoke(null);}catch(Throwable t){return false;} }\n public static boolean isUsingRenderPack(){ if(!iris())return false; try{Class<?> c=Class.forName("net.irisshaders.iris.api.v0.IrisApi"); Object i=c.getMethod("getInstance").invoke(null); return (Boolean)c.getMethod("isShaderPackInUse").invoke(i);}catch(Throwable t){return false;} }\n public static boolean assignScopePipelineToHand(RenderPipeline p,String n){ if(!iris())return false; try{Class<?> a=Class.forName("net.irisshaders.iris.api.v0.IrisApi"); Object i=a.getMethod("getInstance").invoke(null); Class<?> pr=Class.forName("net.irisshaders.iris.api.v0.IrisProgram"); Object hand=Enum.valueOf((Class)pr.asSubclass(Enum.class),"HAND"); a.getMethod("assignPipeline",RenderPipeline.class,pr).invoke(i,p,hand); return true;}catch(Throwable t){return false;} }\n public static void assignCommonEntityPipelinesToHandIfNeeded(){} public static boolean shouldDisableScopeMaskUnderShaderPack(){return false;}\n public static boolean isHandRendererActive(){ if(!isUsingRenderPack())return false; try{Class<?> c=Class.forName("net.irisshaders.iris.pathways.HandRenderer"); Object i=c.getField("INSTANCE").get(null); return (Boolean)c.getMethod("isActive").invoke(i);}catch(Throwable t){return false;} }\n public static boolean endBatch(Object o){return false;} public static boolean endBatch(SubmitNodeCollector c){return false;}\n}\n''')

# Mojang 26.2 renames for retained NeoForge-specific code.
for p in java_root.rglob("*.java"):
    text=p.read_text(encoding="utf-8",errors="ignore")
    text=text.replace("net.minecraft.resources.ResourceLocation","net.minecraft.resources.Identifier")
    text=re.sub(r"\bResourceLocation\b","Identifier",text)
    text=text.replace("import net.minecraft.Util;","import net.minecraft.util.Util;")
    text=text.replace("import net.minecraft.client.renderer.RenderType;","import net.minecraft.client.renderer.rendertype.RenderType;")
    text=text.replace("import net.minecraft.client.renderer.LightTexture;","import net.minecraft.client.renderer.Lightmap;")
    text=re.sub(r"\bLightTexture\b","Lightmap",text)
    p.write_text(text,encoding="utf-8")

# KubeJS optional bridge.
kube=java_root/"com/tacz/guns/api/event/common/KubeJSGunEventPoster.java"
if kube.exists(): kube.write_text('''package com.tacz.guns.api.event.common;\nimport net.neoforged.bus.api.Event; public interface KubeJSGunEventPoster<E extends Event>{ default void postEventToKubeJS(E e){} default void postClientEventToKubeJS(E e){} default void postServerEventToKubeJS(E e){} }\n''',encoding="utf-8")

# Keep only stable public constants/hooks from the old integration registry.
compat=java_root/"com/tacz/guns/init/CompatRegistry.java"
if compat.exists():
    compat.write_text('''package com.tacz.guns.init; public final class CompatRegistry { public static final String IRIS="iris"; private CompatRegistry(){} public static void init(){} public static void initClient(){} }\n''',encoding="utf-8")

p=root/"src/main/resources/kubejs.plugins.txt"
if p.exists(): p.write_text("",encoding="utf-8")
print("Applied TACZ NeoForge 26.2 migration pass 3")
