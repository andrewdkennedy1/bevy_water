use bevy::prelude::*;
use bevy::window::PrimaryWindow;

#[derive(Resource, Default)]
pub struct Pointer {
    pub down: bool,
    pub world: Vec2,
    pub delta: Vec2,
}

pub fn update_pointer(
    buttons: Res<ButtonInput<MouseButton>>,
    windows: Query<&Window, With<PrimaryWindow>>,
    cam_q: Query<(&Camera, &GlobalTransform)>,
    mut pointer: ResMut<Pointer>,
) {
    let window = if let Some(w) = windows.iter().next() { w } else { return };
    let (camera, cam_xform) = if let Some(c) = cam_q.iter().next() { c } else { return };

    let prev = pointer.world;
    pointer.down = buttons.pressed(MouseButton::Left);

    if let Some(cursor) = window.cursor_position() {
        if let Ok(world) = camera.viewport_to_world_2d(cam_xform, cursor) {
            pointer.world = world;
        }
    }
    pointer.delta = pointer.world - prev;
}
