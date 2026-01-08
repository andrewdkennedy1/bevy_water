use bevy::prelude::*;
use std::collections::HashMap;

type FastHashMap<K, V> = HashMap<K, V>;
use crate::input::Pointer;
use crate::constants::BOUNDS_HALF;
use bevy_egui::{egui, EguiContexts};

const PARTICLES_X: usize = 45;
const PARTICLES_Y: usize = 25;
const PARTICLE_SPACING: f32 = 13.0;
const PARTICLE_RADIUS: f32 = 5.0;
const DT: f32 = 1.0 / 60.0;
const BOUNCE: f32 = 0.45;
const DRAG_RADIUS: f32 = 70.0;

#[derive(Resource)]
pub struct SimSettings {
    pub gravity: f32,
    pub viscosity: f32,
    pub stiffness: f32,
    pub rest_density: f32,
    pub h: f32,
    pub drag_strength: f32,
}

impl Default for SimSettings {
    fn default() -> Self {
        Self {
            gravity: -100.0,
            viscosity: 0.15,
            stiffness: 40.0,
            rest_density: 1.0,
            h: 15.0,
            drag_strength: 20.0,
        }
    }
}

#[derive(Component)]
pub struct ParticleVisual;

#[derive(Component, Default)]
struct WaterParticle;

#[derive(Component)]
struct ParticleIndex(usize);

#[derive(Component, Copy, Clone, Default)]
struct Pos(Vec2);

#[derive(Component, Copy, Clone, Default)]
struct Vel(Vec2);

#[derive(Component, Copy, Clone, Default)]
struct Density(f32);

#[derive(Component)]
struct ParticleMat(Handle<ColorMaterial>);

#[derive(Resource, Default)]
struct SpatialHash {
    map: FastHashMap<(i32, i32), Vec<usize>>,
    positions: Vec<Vec2>,
    velocities: Vec<Vec2>,
    densities: Vec<f32>,
    pressures: Vec<f32>,
}

pub struct ParticleFluidPlugin;

impl Plugin for ParticleFluidPlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<SimSettings>()
            .init_resource::<SpatialHash>()
            .add_systems(Startup, particle_setup)
            .add_systems(
                Update,
                (
                    ui_system,
                    build_spatial_hash,
                    compute_density_and_pressure.after(build_spatial_hash),
                    compute_forces.after(compute_density_and_pressure),
                    integrate_and_render.after(compute_forces),
                ),
            );
    }
}

fn ui_system(
    mut contexts: EguiContexts, 
    mut settings: ResMut<SimSettings>, 
    mut q: Query<&mut Vel, With<WaterParticle>>,
    mut frame: Local<u32>
) {
    // Wait 2 frames for Egui fonts to initialize
    *frame += 1;
    if *frame < 3 { return; }
    
    let Ok(ctx) = contexts.ctx_mut() else { return; };
    
    egui::Window::new("Fluid Controls").show(ctx, |ui| {
        ui.add(egui::Slider::new(&mut settings.gravity, -500.0..=500.0).text("Gravity"));
        ui.add(egui::Slider::new(&mut settings.viscosity, 0.0..=1.0).text("Viscosity"));
        ui.add(egui::Slider::new(&mut settings.stiffness, 0.0..=200.0).text("Stiffness"));
        ui.add(egui::Slider::new(&mut settings.rest_density, 0.5..=5.0).text("Rest Density"));
        ui.add(egui::Slider::new(&mut settings.h, 5.0..=40.0).text("Part. Influence (H)"));
        ui.add(egui::Slider::new(&mut settings.drag_strength, 0.0..=100.0).text("Mouse Strength"));
        
        ui.separator();
        ui.horizontal(|ui| {
            if ui.button("Zero G").clicked() { settings.gravity = 0.0; }
            if ui.button("Reset").clicked() { *settings = SimSettings::default(); }
            if ui.button("Freeze").clicked() {
                for mut v in &mut q { v.0 = Vec2::ZERO; }
            }
        });
    });
}

fn particle_setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<ColorMaterial>>,
) {
    let circle = Mesh::from(Circle::new(PARTICLE_RADIUS));
    let circle_handle = meshes.add(circle);

    let origin = Vec2::new(-BOUNDS_HALF.x * 0.65, -BOUNDS_HALF.y * 0.35);
    let mut idx = 0;
    for y in 0..PARTICLES_Y {
        for x in 0..PARTICLES_X {
            let p = origin + Vec2::new(x as f32, y as f32) * PARTICLE_SPACING;
            let mat = ColorMaterial::from(Color::srgb(0.08, 0.25, 0.85));
            let mat_handle = materials.add(mat);

            commands.spawn((
                ParticleVisual,
                WaterParticle,
                ParticleIndex(idx),
                Pos(p),
                Vel(Vec2::ZERO),
                Density(1.0),
                ParticleMat(mat_handle.clone()),
                Mesh2d(circle_handle.clone()),
                MeshMaterial2d(mat_handle),
                Transform::from_translation(p.extend(0.0)),
            ));
            idx += 1;
        }
    }

    commands.spawn((
        ParticleVisual,
        Sprite {
            color: Color::srgba(1.0, 1.0, 1.0, 0.06),
            custom_size: Some(BOUNDS_HALF * 2.0),
            ..default()
        },
        Transform::from_translation(Vec3::new(0.0, 0.0, -1.0)),
    ));
}

#[inline]
fn cell_of(p: Vec2, h: f32) -> (i32, i32) {
    ((p.x / h).floor() as i32, (p.y / h).floor() as i32)
}

fn build_spatial_hash(
    settings: Res<SimSettings>,
    pointer: Res<Pointer>,
    mut contexts: EguiContexts,
    mut sh: ResMut<SpatialHash>,
    q: Query<(&Pos, &Vel, &ParticleIndex), With<WaterParticle>>
) {
    sh.map.clear();
    let n = q.iter().count();
    sh.positions.resize(n, Vec2::ZERO);
    sh.velocities.resize(n, Vec2::ZERO);

    for (p, v, idx) in q.iter() {
        let i = idx.0;
        sh.positions[i] = p.0;
        sh.velocities[i] = v.0;
        let c = cell_of(p.0, settings.h);
        sh.map.entry(c).or_default().push(i);
    }

    // Skip interaction if pointer is over UI
    if contexts.ctx_mut().map_or(false, |ctx| ctx.is_pointer_over_area()) {
        return;
    }

    let cursor = pointer.world;
    let h = settings.h;
    let inv_radius = 1.0 / DRAG_RADIUS;
    let is_moving = pointer.delta.length_squared() > 1.0;
    
    // Use split borrows to avoid conflict between map lookup and velocity update
    let SpatialHash { map, positions, velocities, .. } = &mut *sh;

    if is_moving || pointer.down {
        // Only check nearby cells for interaction
        let (cx, cy) = cell_of(cursor, h);
        let range = (DRAG_RADIUS / h).ceil() as i32;
        
        for dy in -range..=range {
            for dx in -range..=range {
                if let Some(indices) = map.get(&(cx + dx, cy + dy)) {
                    for &i in indices {
                        let d = positions[i] - cursor;
                        let r2 = d.length_squared();
                        if r2 < DRAG_RADIUS * DRAG_RADIUS {
                            let r = r2.sqrt().max(0.001);
                            let factor = (1.0 - r * inv_radius).powf(2.0);
                            if is_moving {
                                velocities[i] += pointer.delta * factor * settings.drag_strength * 5.0;
                            }
                            if pointer.down {
                                 let push = d * (factor * settings.drag_strength * 10.0 * DT / r);
                                 velocities[i] += push;
                            }
                        }
                    }
                }
            }
        }
    }
}

#[inline]
fn poly6(r: f32, h: f32) -> f32 {
    if r >= 0.0 && r <= h {
        let x = 1.0 - (r / h);
        x * x * x
    } else {
        0.0
    }
}

#[inline]
fn spiky_grad(r_vec: Vec2, h: f32) -> Vec2 {
    let r = r_vec.length();
    if r > 0.0001 && r <= h {
        let x = 1.0 - (r / h);
        (-3.0 * x * x / h) * (r_vec / r)
    } else {
        Vec2::ZERO
    }
}

#[inline]
fn viscosity_lap(r: f32, h: f32) -> f32 {
    if r >= 0.0 && r <= h { 1.0 - (r / h) } else { 0.0 }
}

fn neighbors<'a>(sh: &'a SpatialHash, p: Vec2, h: f32) -> impl Iterator<Item = usize> + 'a {
    let (cx, cy) = cell_of(p, h);
    (-1..=1).flat_map(move |dy| {
        (-1..=1).filter_map(move |dx| {
            sh.map.get(&(cx + dx, cy + dy))
        }).flat_map(|v| v.iter().copied())
    })
}

fn compute_density_and_pressure(
    settings: Res<SimSettings>,
    mut sh: ResMut<SpatialHash>,
    mut q: Query<(&mut Density, &ParticleIndex), With<WaterParticle>>
) {
    let n = sh.positions.len();
    sh.densities.resize(n, 0.0);
    sh.pressures.resize(n, 0.0);

    let h = settings.h;
    q.iter_mut().for_each(|(mut d, idx)| {
        let i = idx.0;
        let pi = sh.positions[i];
        let mut dens = 0.0;
        for j in neighbors(&sh, pi, h) {
            let r = (pi - sh.positions[j]).length();
            dens += poly6(r, h);
        }
        d.0 = dens.max(0.0001);
    });

    for (d, idx) in q.iter() {
        let i = idx.0;
        sh.densities[i] = d.0;
        sh.pressures[i] = settings.stiffness * (d.0 - settings.rest_density);
    }
}

fn compute_forces(
    settings: Res<SimSettings>,
    sh: Res<SpatialHash>,
    mut q: Query<(&mut Vel, &Density, &ParticleIndex), With<WaterParticle>>
) {
    let h = settings.h;
    let gravity = Vec2::new(0.0, settings.gravity);

    q.iter_mut().for_each(|(mut v, d, idx)| {
        let i = idx.0;
        let pi = sh.positions[i];
        let vi = sh.velocities[i];
        let di = d.0;
        let pi_pressure = settings.stiffness * (di - settings.rest_density);

        let mut f_pressure = Vec2::ZERO;
        let mut f_visc = Vec2::ZERO;

        for j in neighbors(&sh, pi, h) {
            if i == j { continue; }
            let pj = sh.positions[j];
            let vj = sh.velocities[j];
            let dj = sh.densities[j];
            let pj_pressure = sh.pressures[j];

            let rij = pi - pj;
            let r = rij.length();
            if r > h || r < 0.0001 { continue; }

            let grad = spiky_grad(rij, h);
            let press_term = (pi_pressure + pj_pressure) / (2.0 * dj.max(0.0001));
            f_pressure += -press_term * grad;

            let lap = viscosity_lap(r, h);
            f_visc += settings.viscosity * (vj - vi) * lap;
        }

        let force = f_pressure + f_visc + gravity - 0.9 * vi;
        v.0 += force * DT;
    });
}

fn integrate_and_render(
    mut materials: ResMut<Assets<ColorMaterial>>,
    mut q: Query<(&mut Transform, &mut Pos, &mut Vel, &ParticleMat), With<WaterParticle>>,
) {
    q.iter_mut().for_each(|(mut t, mut p, mut v, mat_comp)| {
        let mut pos = p.0;
        let mut vel = v.0;
        pos += vel * DT;

        if pos.x < -BOUNDS_HALF.x { pos.x = -BOUNDS_HALF.x; vel.x = -vel.x * BOUNCE; }
        if pos.x > BOUNDS_HALF.x { pos.x = BOUNDS_HALF.x; vel.x = -vel.x * BOUNCE; }
        if pos.y < -BOUNDS_HALF.y { pos.y = -BOUNDS_HALF.y; vel.y = -vel.y * BOUNCE; }
        if pos.y > BOUNDS_HALF.y { pos.y = BOUNDS_HALF.y; vel.y = -vel.y * BOUNCE; }

        p.0 = pos;
        v.0 = vel;
        t.translation.x = pos.x;
        t.translation.y = pos.y;

        // Update visuals in the same loop
        if let Some(mat) = materials.get_mut(&mat_comp.0) {
            let speed = vel.length();
            let foam = (speed / 280.0).clamp(0.0, 1.0);
            let energy = (speed / 450.0).clamp(0.0, 1.0);
            
            mat.color = Color::srgb(
                0.06 + 0.35 * foam + 0.60 * energy * energy,
                0.22 + 0.65 * foam + 0.15 * energy,
                0.75 + 0.30 * foam - 0.20 * energy,
            );
        }
    });
}
