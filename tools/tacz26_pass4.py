#!/usr/bin/env python3
from pathlib import Path
import json
import re
import sys

root = Path(sys.argv[1]).resolve()
fabric = root.parent / "fabric26"
java = root / "src/main/java"
res = root / "src/main/resources"
fabjava = fabric / "src/main/java"


def copy_fabric(rel: str) -> Path:
    src = fabjava / rel
    dst = java / rel
    if not src.exists():
        raise FileNotFoundError(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    text = src.read_text(encoding="utf-8", errors="ignore")
    text = text.replace("import net.fabricmc.api.EnvType;\n", "")
    text = text.replace("import net.fabricmc.api.Environment;\n", "")
    text = re.sub(r"\s*@Environment\(EnvType\.(?:CLIENT|SERVER)\)\s*", "\n", text)
    dst.write_text(text, encoding="utf-8")
    return dst


def replace_file(path: Path, pairs):
    text = path.read_text(encoding="utf-8", errors="ignore")
    for a, b in pairs:
        text = text.replace(a, b)
    path.write_text(text, encoding="utf-8")


# ---------------------------------------------------------------------------
# 26.2 GUI migration. The 26.2 Fabric source already contains the vanilla GUI
# rewrite. Only the transport and Screen accessor are loader-specific.
# ---------------------------------------------------------------------------
for rel in [
    "com/tacz/guns/client/gui/GunSmithTableScreen.java",
    "com/tacz/guns/client/gui/GunRefitScreen.java",
]:
    p = copy_fabric(rel)
    t = p.read_text(encoding="utf-8")
    t = t.replace("import cn.sh1rocu.tacz.mixin.accessor.ScreenAccessor;",
                  "import com.tacz.guns.mixin.client.ScreenAccessor;")
    t = t.replace("import net.fabricmc.fabric.api.client.networking.v1.ClientPlayNetworking;",
                  "import net.neoforged.neoforge.network.PacketDistributor;")
    t = t.replace("ClientPlayNetworking.send(", "PacketDistributor.sendToServer(")
    p.write_text(t, encoding="utf-8")

# Local equivalent of the Fabric accessor, registered in the existing TACZ
# Mixin config so no Fabric helper package is needed at runtime.
accessor = java / "com/tacz/guns/mixin/client/ScreenAccessor.java"
accessor.parent.mkdir(parents=True, exist_ok=True)
accessor.write_text('''package com.tacz.guns.mixin.client;\n\nimport net.minecraft.client.gui.components.Renderable;\nimport net.minecraft.client.gui.screens.Screen;\nimport org.spongepowered.asm.mixin.Mixin;\nimport org.spongepowered.asm.mixin.gen.Accessor;\nimport java.util.List;\n\n@Mixin(Screen.class)\npublic interface ScreenAccessor {\n    @Accessor("renderables")\n    List<Renderable> tacz$getRenderables();\n}\n''', encoding="utf-8")

mixins = res / "tacz.mixins.json"
if mixins.exists():
    data = json.loads(mixins.read_text(encoding="utf-8"))
    data["compatibilityLevel"] = "JAVA_25"
    client = data.setdefault("client", [])
    if "client.ScreenAccessor" not in client:
        client.append("client.ScreenAccessor")
    mixins.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

# ---------------------------------------------------------------------------
# Client configs added by the 26.2 port. NeoForge uses ModConfigSpec directly.
# ---------------------------------------------------------------------------
for rel in [
    "com/tacz/guns/config/client/SoundConfig.java",
    "com/tacz/guns/config/client/ResourceConfig.java",
]:
    p = copy_fabric(rel)
    replace_file(p, [
        ("net.minecraftforge.common.ForgeConfigSpec", "net.neoforged.neoforge.common.ModConfigSpec"),
        ("ForgeConfigSpec", "ModConfigSpec"),
    ])

client_config = java / "com/tacz/guns/config/ClientConfig.java"
if client_config.exists():
    text = client_config.read_text(encoding="utf-8")
    if "ResourceConfig.init(builder);" not in text:
        text = text.replace("RenderConfig.init(builder);",
                            "RenderConfig.init(builder);\n        ResourceConfig.init(builder);\n        SoundConfig.init(builder);")
    if "import com.tacz.guns.config.client.*;" not in text:
        text = text.replace("import com.tacz.guns.config.client.KeyConfig;\nimport com.tacz.guns.config.client.RenderConfig;\nimport com.tacz.guns.config.client.ZoomConfig;",
                            "import com.tacz.guns.config.client.*;")
    client_config.write_text(text, encoding="utf-8")

sound_mgr = java / "com/tacz/guns/client/sound/SoundPlayManager.java"
if sound_mgr.exists():
    text = sound_mgr.read_text(encoding="utf-8")
    text = re.sub(r"ModSounds\.GUN(?!\.get\(\))", "ModSounds.GUN.get()", text)
    sound_mgr.write_text(text, encoding="utf-8")

# ---------------------------------------------------------------------------
# Minecart target: take Mojang 26.2 behavior from the working Fabric port, then
# restore NeoForge registrations and event bus semantics.
# ---------------------------------------------------------------------------
target = copy_fabric("com/tacz/guns/entity/TargetMinecart.java")
text = target.read_text(encoding="utf-8")
text = text.replace("import cn.sh1rocu.tacz.api.LogicalSide;", "import net.neoforged.fml.LogicalSide;\nimport net.neoforged.neoforge.common.NeoForge;")
text = text.replace("import cn.sh1rocu.tacz.api.extension.IMinecart;\n", "")
text = text.replace("implements ITargetEntity, IMinecart", "implements ITargetEntity")
text = text.replace("ModSounds.TARGET_HIT,", "ModSounds.TARGET_HIT.get(),")
text = text.replace("new ItemStack(ModItems.TARGET_MINECART)", "new ItemStack(ModItems.TARGET_MINECART.get())")
text = text.replace("return ModItems.TARGET_MINECART;", "return ModItems.TARGET_MINECART.get();")
text = text.replace("return ModBlocks.TARGET.defaultBlockState();", "return ModBlocks.TARGET.get().defaultBlockState();")
text = text.replace("EntityHurtByGunEvent.Post event = new EntityHurtByGunEvent.Post(projectile, this, player, projectile.getGunId(), projectile.getGunDisplayId(), damage, Pair.of(source, source), isHeadshot, headshotMultiplier, LogicalSide.SERVER);\n                EntityHurtByGunEvent.POST.invoker().post(event);",
                    "EntityHurtByGunEvent.Post event = new EntityHurtByGunEvent.Post(projectile, this, player, projectile.getGunId(), projectile.getGunDisplayId(), damage, Pair.of(source, source), isHeadshot, headshotMultiplier, LogicalSide.SERVER);\n                NeoForge.EVENT_BUS.post(event);")
# Fabric's IMinecart mixin supplies this extension point. NeoForge does not need
# the shim for compilation; remove the method rather than pretending it overrides
# a vanilla API.
text = re.sub(r"\n\s*@Override\n\s*public boolean tacz\$canBeRidden\(\) \{\n\s*return false;\n\s*\}\n", "\n", text)
target.write_text(text, encoding="utf-8")

# ---------------------------------------------------------------------------
# Reconstruct 1.21's interpolatable walk distance on 26.2. This is a tiny,
# loader-neutral interface + NeoForge mixin, copied from the proven Fabric port.
# ---------------------------------------------------------------------------
copy_fabric("cn/sh1rocu/tacz/api/extension/IMoveDistTracker.java")
ctx = copy_fabric("com/tacz/guns/client/animation/statemachine/GunAnimationStateContext.java")
move_mixin = java / "com/tacz/guns/mixin/common/EntityMoveDistMixin.java"
move_mixin.parent.mkdir(parents=True, exist_ok=True)
move_mixin.write_text('''package com.tacz.guns.mixin.common;\n\nimport cn.sh1rocu.tacz.api.extension.IMoveDistTracker;\nimport net.minecraft.world.entity.Entity;\nimport org.spongepowered.asm.mixin.Mixin;\nimport org.spongepowered.asm.mixin.Unique;\nimport org.spongepowered.asm.mixin.injection.At;\nimport org.spongepowered.asm.mixin.injection.Inject;\nimport org.spongepowered.asm.mixin.injection.callback.CallbackInfo;\n\n@Mixin(Entity.class)\npublic abstract class EntityMoveDistMixin implements IMoveDistTracker {\n    @Unique private float tacz$moveDistO;\n    @Unique private boolean tacz$moveDistInit;\n    @Unique @Override public float tacz$getMoveDistO() {\n        return this.tacz$moveDistInit ? this.tacz$moveDistO : ((Entity)(Object)this).moveDist;\n    }\n    @Inject(method = "tick", at = @At("HEAD"))\n    private void tacz$captureMoveDistO(CallbackInfo ci) {\n        this.tacz$moveDistO = ((Entity)(Object)this).moveDist;\n        this.tacz$moveDistInit = true;\n    }\n}\n''', encoding="utf-8")
if mixins.exists():
    data = json.loads(mixins.read_text(encoding="utf-8"))
    common = data.setdefault("mixins", [])
    if "common.EntityMoveDistMixin" not in common:
        common.append("common.EntityMoveDistMixin")
    mixins.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

# ---------------------------------------------------------------------------
# 26.2 item rendering no longer uses the old BEWLR/MultiBufferSource path. Use
# the tested 26.2 ItemModel/SpecialModelRenderer bridge from the working port.
# The registry helper itself is loader-neutral despite its historical package.
# ---------------------------------------------------------------------------
for rel in [
    "cn/sh1rocu/tacz/compat/fabric/BuiltinItemRendererRegistry.java",
    "com/tacz/guns/client/renderer/item/TaczDynamicItemModel.java",
    "com/tacz/guns/client/renderer/item/AttachmentItemRenderer.java",
    "com/tacz/guns/client/renderer/item/AmmoItemRenderer.java",
    "com/tacz/guns/client/renderer/item/GunSmithTableItemRenderer.java",
]:
    copy_fabric(rel)

# PlayerAnimator is optional in the core build, but ClientSetupEvent passes a
# method reference. Preserve the exact Consumer target type so javac can infer it.
pa = java / "com/tacz/guns/compat/playeranimator/PlayerAnimatorCompat.java"
if pa.exists():
    text = pa.read_text(encoding="utf-8")
    text = text.replace("public static void registerReloadListener(Object o){}",
                        "public static void registerReloadListener(java.util.function.Consumer<net.minecraft.server.packs.resources.PreparableReloadListener> register){}")
    pa.write_text(text, encoding="utf-8")

print("Applied TACZ NeoForge 26.2 targeted migration pass 4")
