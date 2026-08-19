#!/usr/bin/env python3
from pathlib import Path
import re
import sys

root = Path(sys.argv[1]).resolve()
fabric = root.parent / "fabric26"
java = root / "src/main/java"
fabjava = fabric / "src/main/java"


def clean_env(text: str) -> str:
    text = text.replace("import net.fabricmc.api.EnvType;\n", "")
    text = text.replace("import net.fabricmc.api.Environment;\n", "")
    return re.sub(r"\s*@Environment\(EnvType\.(?:CLIENT|SERVER)\)\s*", "\n", text)


def copy_rel(rel: str) -> Path:
    src=fabjava/rel; dst=java/rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(clean_env(src.read_text(encoding='utf-8',errors='ignore')),encoding='utf-8')
    return dst

# Pass 3 was intentionally conservative and retained a 1.21.1 file whenever the
# NeoForge version itself used NeoForge APIs. For 26.2 that leaves many files on
# obsolete Mojang rendering/world APIs. If the 26.2 reference version is truly
# loader-neutral, it is safe to prefer it regardless of what the old file used.
common_neo={p.relative_to(java).as_posix():p for p in java.rglob('*.java')}
common_fab={p.relative_to(fabjava).as_posix():p for p in fabjava.rglob('*.java')}
copied=0
for rel in sorted(set(common_neo)&set(common_fab)):
    text=common_fab[rel].read_text(encoding='utf-8',errors='ignore')
    bad=False
    for line in text.splitlines():
        s=line.strip()
        if not s.startswith('import '): continue
        imp=s.removeprefix('import ').removesuffix(';')
        if imp.startswith(('net.fabricmc.fabric.api','net.fabricmc.loader.api','cn.sh1rocu.')):
            bad=True; break
    if not bad:
        common_neo[rel].write_text(clean_env(text),encoding='utf-8')
        copied+=1
print(f'Pass5: refreshed {copied} loader-neutral files from the MC26.2 source')

# 26.2 inventory abstraction used by the working port is local, small and has no
# Fabric runtime dependency. Bring it over so gun ammo inventory logic no longer
# depends on the removed 1.21 NeoForge capability surface.
for src in (fabjava/'cn/sh1rocu/tacz/util/itemhandler').rglob('*.java'):
    rel=src.relative_to(fabjava).as_posix()
    copy_rel(rel)
copy_rel('cn/sh1rocu/tacz/api/extension/IItem.java')

# Dynamic-renderer item implementations. Their only loader-looking imports are
# the local registry and IItem bridge already provided by pass4/pass5.
for rel in [
    'com/tacz/guns/api/item/gun/AbstractGunItem.java',
    'com/tacz/guns/item/AttachmentItem.java',
    'com/tacz/guns/item/AmmoItem.java',
    'com/tacz/guns/item/GunSmithTableItem.java',
]:
    copy_rel(rel)

# ---------------------------------------------------------------------------
# Kinetic bullet: use the validated Mojang 26.2 implementation, retaining
# NeoForge's native complex-spawn packet and event bus. NeoForge still supplies
# Entity#getPersistentData(), so the temporary Fabric entity-data shim is not
# needed here.
# ---------------------------------------------------------------------------
p=copy_rel('com/tacz/guns/entity/EntityKineticBullet.java')
t=p.read_text(encoding='utf-8')
t=t.replace('import cn.sh1rocu.tacz.api.LogicalSide;', 'import net.neoforged.fml.LogicalSide;')
t=t.replace('import cn.sh1rocu.tacz.api.extension.IEntityAdditionalSpawnData;\n', 'import net.neoforged.neoforge.entity.IEntityWithComplexSpawn;\n')
t=t.replace('import cn.sh1rocu.tacz.api.extension.IEntityPersistentData;\n', 'import net.neoforged.neoforge.common.NeoForge;\n')
t=t.replace('import net.minecraft.network.FriendlyByteBuf;', 'import net.minecraft.network.RegistryFriendlyByteBuf;')
t=t.replace('import net.minecraft.network.protocol.Packet;\n', '')
t=t.replace('import net.minecraft.network.protocol.game.ClientGamePacketListener;\n', '')
t=t.replace('import net.minecraft.server.level.ServerEntity;\n', '')
t=t.replace('implements IEntityAdditionalSpawnData', 'implements IEntityWithComplexSpawn')
# NeoForge's spawn extension generates the advanced spawn payload itself.
t=re.sub(r'\n\s*@Override\n\s*public @NotNull Packet<ClientGamePacketListener> getAddEntityPacket\(ServerEntity entity\) \{\n\s*return IEntityAdditionalSpawnData\.getEntitySpawningPacket\(this\);\n\s*\}\n', '\n', t)
t=t.replace('writeSpawnData(FriendlyByteBuf buffer)', 'writeSpawnData(RegistryFriendlyByteBuf buffer)')
t=t.replace('readSpawnData(FriendlyByteBuf additionalData)', 'readSpawnData(RegistryFriendlyByteBuf additionalData)')
t=t.replace('((IEntityPersistentData) this).tacz$getPersistentData()', 'this.getPersistentData()')
# Restore the cancellable NeoForge event semantics from the native port.
t=t.replace('EntityHurtByGunEvent.PRE.invoker().post(preEvent);', 'if (NeoForge.EVENT_BUS.post(preEvent).isCanceled()) { return; }')
t=t.replace('EntityKillByGunEvent.CALLBACK.invoker().post(killByGunEvent);', 'NeoForge.EVENT_BUS.post(killByGunEvent);')
t=t.replace('EntityHurtByGunEvent.POST.invoker().post(hurtByGunEvent);', 'NeoForge.EVENT_BUS.post(hurtByGunEvent);')
t=t.replace('AmmoHitBlockEvent.CALLBACK.invoker().post(ammoHitBlockEvent);', 'if (NeoForge.EVENT_BUS.post(ammoHitBlockEvent).isCanceled()) { return; }')
p.write_text(t,encoding='utf-8')

# ---------------------------------------------------------------------------
# 26.2 NeoForge DeferredRegister factories now inject registry-aware properties.
# Use registerBlock/registerItem for classes whose constructors take Properties.
# ---------------------------------------------------------------------------
blocks=java/'com/tacz/guns/init/ModBlocks.java'
if blocks.exists():
    t=blocks.read_text(encoding='utf-8')
    for name,ctor in [
        ('gun_smith_table','GunSmithTableBlockB'),('workbench_a','GunSmithTableBlockA'),
        ('workbench_b','GunSmithTableBlockB'),('workbench_c','GunSmithTableBlockC'),
        ('target','TargetBlock'),('statue','StatueBlock')]:
        t=t.replace(f'BLOCKS.register("{name}", {ctor}::new)', f'BLOCKS.registerBlock("{name}", {ctor}::new)')
    blocks.write_text(t,encoding='utf-8')

items=java/'com/tacz/guns/init/ModItems.java'
if items.exists():
    t=items.read_text(encoding='utf-8')
    t=t.replace('ITEMS.register("modern_kinetic_gun", ModernKineticGunItem::new)', 'ITEMS.registerItem("modern_kinetic_gun", ModernKineticGunItem::new)')
    t=t.replace('ITEMS.register("target_minecart", TargetMinecartItem::new)', 'ITEMS.registerItem("target_minecart", TargetMinecartItem::new)')
    t=t.replace('ITEMS.register("gun_smith_table", () -> new DefaultTableItem(ModBlocks.GUN_SMITH_TABLE.get()))', 'ITEMS.registerItem("gun_smith_table", props -> new DefaultTableItem(ModBlocks.GUN_SMITH_TABLE.get(), props))')
    items.write_text(t,encoding='utf-8')

# 26.2 EventBusSubscriber no longer selects a separate MOD bus in the annotation.
for p in java.rglob('*.java'):
    t=p.read_text(encoding='utf-8',errors='ignore')
    t=t.replace('@EventBusSubscriber(bus = EventBusSubscriber.Bus.MOD, ', '@EventBusSubscriber(')
    t=t.replace('@EventBusSubscriber(bus = EventBusSubscriber.Bus.MOD)', '@EventBusSubscriber')
    p.write_text(t,encoding='utf-8')

print('Applied TACZ NeoForge 26.2 migration pass 5')
