mod input;
mod constants;
mod particles;

use bevy::prelude::*;
use bevy_egui::EguiPlugin;
use crate::input::{Pointer, update_pointer};
use crate::particles::ParticleFluidPlugin;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_plugins(EguiPlugin::default())
        .init_resource::<Pointer>()
        .add_systems(Startup, setup_camera)
        .add_systems(Update, update_pointer)
        .add_plugins(ParticleFluidPlugin)
        .run();
}

fn setup_camera(mut commands: Commands) {
    commands.spawn(Camera2d::default());
}
