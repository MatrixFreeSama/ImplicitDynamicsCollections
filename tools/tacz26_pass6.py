#!/usr/bin/env python3
from pathlib import Path
import json
import re
import shutil
import subprocess
import sys

root = Path(sys.argv[1]).resolve()
fabric = root.parent / "fabric26"
java = root / "src/main/java"
fabjava = fabric / "src/main/java"
res = root / "src/main/resources"


def git_restore(rel: str):
    subprocess.run(["git", "-C", str(root), "checkout", "HEAD", "--", rel], check=True)


def clean_env(text: str) -> str:
    text = text.replace("import net.fabricmc.api.EnvType;\n", "")
    text = text.replace("import net.fabricmc.api.Environment;\n", "")
    return re.sub(r"\s*@Environment\(EnvType\.(?:CLIENT|SERVER)\)\s*", "\n", text)


def copy_fabric(rel: str) -> Path:
    src = fabjava / rel
    dst = java / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(clean_env(src.read_text(encoding="utf-8", errors="ignore")), encoding="utf-8")
    return dst


def modernize(text: str) -> str:
    text = text.replace("net.minecraft.resources.ResourceLocation", "net.minecraft.resources.Identifier")
    text = re.sub(r"\bResourceLocation\b", "Identifier", text)
    text = text.replace("import net.minecraft.Util;", "import net.minecraft.util.Util;")
    text = text.replace("import net.minecraft.client.renderer.RenderType;", "import net.minecraft.client.renderer.rendertype.RenderType;")
    text = text.replace("import net.minecraft.client.renderer.LightTexture;", "import net.minecraft.client.renderer.Lightmap;")
    text = re.sub(r"\bLightTexture\b", "Lightmap", text)
    text = re.sub(r"@EventBusSubscriber\(([^)]*?)\s*,?\s*bus\s*=\s*EventBusSubscriber\.Bus\.MOD\s*,?\s*([^)]*?)\)", lambda m: "@EventBusSubscriber(" + ", ".join(x.strip(" ,") for x in [m.group(1), m.group(2)] if x.strip(" ,")) + ")", text)
    text = text.replace("@EventBusSubscriber(bus = EventBusSubscriber.Bus.MOD)", "@EventBusSubscriber")
    return text

# Pass5 deliberately tried a broad 26.2 source transplant. A number of loader
# boundary files have identical-looking Java but completely different lifecycle
# semantics. Reset those boundaries to the known NeoForge port, then migrate
# them explicitly below.
for rel in [
    "src/main/java/com/tacz/guns/init",
    "src/main/java/com/tacz/guns/client/init",
    "src/main/java/com/tacz/guns/network",
    "src/main/java/com/tacz/guns/config",
    "src/main/java/com/tacz/guns/api/item/gun/AbstractGunItem.java",
    "src/main/java/com/tacz/guns/item/AttachmentItem.java",
    "src/main/java/com/tacz/guns/item/AmmoItem.java",
    "src/main/java/com/tacz/guns/item/GunSmithTableItem.java",
]:
    git_restore(rel)

# Remove the experimental Fabric-side inventory abstraction brought in by pass5.
# The 1.21 NeoForge gun/item core already compiled much further without it.
for rel in [
    "cn/sh1rocu/tacz/util/itemhandler",
    "cn/sh1rocu/tacz/api/extension/IItem.java",
]:
    p = java / rel
    if p.is_dir(): shutil.rmtree(p)
    elif p.exists(): p.unlink()

# Modernize the restored loader-boundary files for Mojang/NeoForge 26.2 names.
for rel in ["com/tacz/guns/init", "com/tacz/guns/client/init", "com/tacz/guns/network", "com/tacz/guns/config"]:
    for p in (java / rel).rglob("*.java"):
        p.write_text(modernize(p.read_text(encoding="utf-8", errors="ignore")), encoding="utf-8")

# ---------------------------------------------------------------------------
# Native NeoForge deferred registration. In 26.2 Block/Item constructors receive
# registry-aware Properties, so use the specialized registerBlock/registerItem
# factories instead of handing those constructors an Identifier.
# ---------------------------------------------------------------------------
blocks = java / "com/tacz/guns/init/ModBlocks.java"
t = blocks.read_text(encoding="utf-8")
for name, ctor in [
    ("gun_smith_table", "GunSmithTableBlockB"), ("workbench_a", "GunSmithTableBlockA"),
    ("workbench_b", "GunSmithTableBlockB"), ("workbench_c", "GunSmithTableBlockC"),
    ("target", "TargetBlock"), ("statue", "StatueBlock")
]:
    t = t.replace(f'BLOCKS.register("{name}", {ctor}::new)', f'BLOCKS.registerBlock("{name}", {ctor}::new)')
blocks.write_text(t, encoding="utf-8")

items = java / "com/tacz/guns/init/ModItems.java"
t = items.read_text(encoding="utf-8")
t = t.replace('ITEMS.register("modern_kinetic_gun", ModernKineticGunItem::new)', 'ITEMS.registerItem("modern_kinetic_gun", ModernKineticGunItem::new)')
t = t.replace('ITEMS.register("gun_smith_table", () -> new DefaultTableItem(ModBlocks.GUN_SMITH_TABLE.get()))', 'ITEMS.registerItem("gun_smith_table", props -> new DefaultTableItem(ModBlocks.GUN_SMITH_TABLE.get(), props))')
t = t.replace('ITEMS.register("target_minecart", TargetMinecartItem::new)', 'ITEMS.registerItem("target_minecart", TargetMinecartItem::new)')
items.write_text(t, encoding="utf-8")

# Optional integrations must never stop the core port compiling.
compat = java / "com/tacz/guns/init/CompatRegistry.java"
compat.write_text('''package com.tacz.guns.init;\n\nimport net.neoforged.fml.ModList;\n\npublic final class CompatRegistry {\n    public static final String IRIS = "iris";\n    public static final String CLOTH_CONFIG = "cloth_config";\n    private CompatRegistry() {}\n    public static void init() {}\n    public static void initClient() {}\n    public static void checkModLoad(String id, Runnable action) { if (ModList.get().isLoaded(id)) action.run(); }\n}\n''', encoding="utf-8")
menu = java / "com/tacz/guns/compat/cloth/MenuIntegration.java"
menu.parent.mkdir(parents=True, exist_ok=True)
menu.write_text('''package com.tacz.guns.compat.cloth;\nimport net.minecraft.client.gui.screens.Screen;\npublic final class MenuIntegration { private MenuIntegration(){} public static Screen getConfigScreen(Screen parent){ return null; } }\n''', encoding="utf-8")

# Newer 26.2 client code expects these two config groups. They are ordinary
# NeoForge ModConfigSpec entries and are safe to add to the existing client spec.
for name in ["SoundConfig", "ResourceConfig"]:
    src = fabjava / f"com/tacz/guns/config/client/{name}.java"
    dst = java / f"com/tacz/guns/config/client/{name}.java"
    text = src.read_text(encoding="utf-8").replace("net.minecraftforge.common.ForgeConfigSpec", "net.neoforged.neoforge.common.ModConfigSpec").replace("ForgeConfigSpec", "ModConfigSpec")
    dst.write_text(text, encoding="utf-8")
client_cfg = java / "com/tacz/guns/config/ClientConfig.java"
t = client_cfg.read_text(encoding="utf-8")
t = t.replace("import com.tacz.guns.config.client.KeyConfig;\nimport com.tacz.guns.config.client.RenderConfig;\nimport com.tacz.guns.config.client.ZoomConfig;", "import com.tacz.guns.config.client.*;")
if "SoundConfig.init(builder);" not in t:
    t = t.replace("RenderConfig.init(builder);", "RenderConfig.init(builder);\n        ResourceConfig.init(builder);\n        SoundConfig.init(builder);")
client_cfg.write_text(t, encoding="utf-8")

# ---------------------------------------------------------------------------
# 26.2 key mapping API. Use the validated 26.2 key implementations, but dispatch
# them from NeoForge's native InputEvent and ClientTickEvent.
# ---------------------------------------------------------------------------
for name in ["AimKey", "ConfigKey", "CrawlKey", "FireSelectKey", "InspectKey", "InteractKey", "MeleeKey", "RefitKey", "ReloadKey", "ShootKey", "TaCZKeyCategory", "ZoomKey"]:
    p = copy_fabric(f"com/tacz/guns/client/input/{name}.java")
    text = p.read_text(encoding="utf-8")
    text = text.replace("import cn.sh1rocu.tacz.api.event.InputEvent;", "import net.neoforged.neoforge.client.event.InputEvent;")
    text = text.replace("import cn.sh1rocu.tacz.api.event.PlayerTickEvent;\n", "")
    text = text.replace("import net.fabricmc.loader.api.FabricLoader;", "import net.neoforged.fml.ModList;")
    text = text.replace("FabricLoader.getInstance().isModLoaded", "ModList.get().isLoaded")
    text = text.replace("import net.fabricmc.fabric.api.client.networking.v1.ClientPlayNetworking;", "import net.neoforged.neoforge.network.PacketDistributor;")
    text = text.replace("ClientPlayNetworking.send(", "PacketDistributor.sendToServer(")
    if name == "ReloadKey":
        text = text.replace("public static void autoReload(PlayerTickEvent.Pre event) {\n        if (!event.getEntity().level().isClientSide())\n            return;\n", "public static void autoReload() {\n")
    p.write_text(text, encoding="utf-8")

subscriber = java / "com/tacz/guns/client/input/ClientInputEvents.java"
subscriber.write_text('''package com.tacz.guns.client.input;\n\nimport com.tacz.guns.GunMod;\nimport net.minecraft.client.Minecraft;\nimport net.neoforged.api.distmarker.Dist;\nimport net.neoforged.bus.api.SubscribeEvent;\nimport net.neoforged.fml.common.EventBusSubscriber;\nimport net.neoforged.neoforge.client.event.ClientTickEvent;\nimport net.neoforged.neoforge.client.event.InputEvent;\n\n@EventBusSubscriber(value = Dist.CLIENT, modid = GunMod.MOD_ID)\npublic final class ClientInputEvents {\n    private ClientInputEvents() {}\n    @SubscribeEvent public static void key(InputEvent.Key e) {\n        ConfigKey.onOpenConfig(e); CrawlKey.onCrawlPress(e); FireSelectKey.onFireSelectPress(e);\n        InspectKey.onInspectPress(e); InteractKey.onInteractPress(e); MeleeKey.onMeleePress(e);\n        RefitKey.onRefitPress(e); ReloadKey.onReloadPress(e); ZoomKey.onZoomKeyPress(e);\n    }\n    @SubscribeEvent public static void mouse(InputEvent.MouseButton.Post e) {\n        AimKey.onAimPress(e); ZoomKey.onZoomMousePress(e);\n    }\n    @SubscribeEvent public static void tick(ClientTickEvent.Post e) {\n        Minecraft mc = Minecraft.getInstance();\n        AimKey.onAimHoldingPreInput(mc); AimKey.cancelAim(mc); ReloadKey.autoReload(); ShootKey.autoShoot(mc, true);\n    }\n}\n''', encoding="utf-8")

# ---------------------------------------------------------------------------
# Crosshair/HUD extraction was rewritten by Mojang in 26.2. Use the proven 26.2
# implementation and remove its Fabric-only render-tick event. Screen state is
# sampled immediately before extracting the layer instead.
# ---------------------------------------------------------------------------
cross = copy_fabric("com/tacz/guns/client/event/RenderCrosshairEvent.java")
t = cross.read_text(encoding="utf-8")
t = t.replace("import cn.sh1rocu.simplebedrockmodel.api.event.RenderTickEvent;\n", "")
t = re.sub(r"\n\s*public static void onRenderTick\(RenderTickEvent event\) \{.*?\n\s*\}\n", "\n", t, flags=re.S)
t = t.replace("LocalPlayer player = Minecraft.getInstance().player;", "isRefitScreen = Minecraft.getInstance().gui.screen() instanceof GunRefitScreen;\n        LocalPlayer player = Minecraft.getInstance().player;", 1)
cross.write_text(t, encoding="utf-8")

# Re-run annotation/name cleanup after all restores/copies.
for p in java.rglob("*.java"):
    p.write_text(modernize(p.read_text(encoding="utf-8", errors="ignore")), encoding="utf-8")

print("Applied TACZ NeoForge 26.2 migration pass 6")
