#!/usr/bin/env python3
from pathlib import Path
import sys

root = Path(sys.argv[1]).resolve()
java = root / "src/main/java"

subscriber = java / "com/tacz/guns/client/input/ClientInputEvents.java"
subscriber.write_text('''package com.tacz.guns.client.input;\n\nimport com.tacz.guns.GunMod;\nimport net.minecraft.client.Minecraft;\nimport net.neoforged.api.distmarker.Dist;\nimport net.neoforged.bus.api.SubscribeEvent;\nimport net.neoforged.fml.common.EventBusSubscriber;\nimport net.neoforged.neoforge.client.event.ClientTickEvent;\nimport net.neoforged.neoforge.client.event.InputEvent;\n\n@EventBusSubscriber(value = Dist.CLIENT, modid = GunMod.MOD_ID)\npublic final class ClientInputEvents {\n    private ClientInputEvents() {}\n\n    @SubscribeEvent\n    public static void key(InputEvent.Key e) {\n        ConfigKey.onOpenConfig(e);\n        CrawlKey.onCrawlPress(e);\n        FireSelectKey.onFireSelectKeyPress(e);\n        InspectKey.onInspectPress(e);\n        InteractKey.onInteractKeyPress(e);\n        MeleeKey.onMeleeKeyPress(e);\n        RefitKey.onRefitPress(e);\n        ReloadKey.onReloadPress(e);\n        ZoomKey.onZoomKeyPress(e);\n    }\n\n    @SubscribeEvent\n    public static void mouse(InputEvent.MouseButton.Post e) {\n        AimKey.onAimPress(e);\n        FireSelectKey.onFireSelectMousePress(e);\n        InteractKey.onInteractMousePress(e);\n        MeleeKey.onMeleeMousePress(e);\n        ZoomKey.onZoomMousePress(e);\n    }\n\n    @SubscribeEvent\n    public static void tick(ClientTickEvent.Post e) {\n        Minecraft mc = Minecraft.getInstance();\n        AimKey.onAimHoldingPreInput(mc);\n        AimKey.cancelAim(mc);\n        ReloadKey.autoReload();\n        ShootKey.autoShoot(mc, true);\n    }\n}\n''', encoding="utf-8")

# FancyModLoader 26.x exposes the process side through FMLEnvironment rather
# than the old static FMLLoader#getDist helper.
gunmod = java / "com/tacz/guns/GunMod.java"
if gunmod.exists():
    text = gunmod.read_text(encoding="utf-8")
    text = text.replace("import net.neoforged.fml.loading.FMLLoader;", "import net.neoforged.fml.loading.FMLEnvironment;")
    text = text.replace("FMLLoader.getDist()", "FMLEnvironment.getDist()")
    gunmod.write_text(text, encoding="utf-8")

print("Applied TACZ NeoForge 26.2 pass 6b")
