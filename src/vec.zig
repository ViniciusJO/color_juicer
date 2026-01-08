const std = @import("std");

pub const Vec2 = [2]u8;
pub const Vec3 = [3]u8;
pub const Vec4 = [4]u8;

pub fn vec3_dist_sz(v1: Vec3, v2: [4]usize) f64 {
    const dx: f64 = @floatFromInt(v1[0] - v2[0]);
    const dy: f64 = @floatFromInt(v1[1] - v2[1]);
    const dz: f64 = @floatFromInt(v1[2] - v2[2]);
    return @sqrt(dx * dx + dy * dy + dz * dz);
}

pub fn vec3_dist(v1: Vec3, v2: Vec3) f64 {
    const dx: f64 = @floatFromInt(@abs(@as(i16, v1[0]) - @as(i16, v2[0])));
    const dy: f64 = @floatFromInt(@abs(@as(i16, v1[1]) - @as(i16, v2[1])));
    const dz: f64 = @floatFromInt(@abs(@as(i16, v1[2]) - @as(i16, v2[2])));
    return @sqrt(dx * dx + dy * dy + dz * dz);
}

pub fn vec4_dist_sz(v1: Vec4, v2: [4]usize) f64 {
    const dx: f64 = @floatFromInt(v1[0] - v2[0]);
    const dy: f64 = @floatFromInt(v1[1] - v2[1]);
    const dz: f64 = @floatFromInt(v1[2] - v2[2]);
    const dw: f64 = @floatFromInt(v1[3] - v2[3]);
    return @sqrt(dx * dx + dy * dy + dz * dz + dw * dw);
}

pub fn vec4_dist(v1: Vec4, v2: Vec4) f64 {
    const dx: f64 = @floatFromInt(@abs(@as(i16, v1[0]) - @as(i16, v2[0])));
    const dy: f64 = @floatFromInt(@abs(@as(i16, v1[1]) - @as(i16, v2[1])));
    const dz: f64 = @floatFromInt(@abs(@as(i16, v1[2]) - @as(i16, v2[2])));
    const dw: f64 = @floatFromInt(@abs(@as(i16, v1[3]) - @as(i16, v2[3])));
    return @sqrt(dx * dx + dy * dy + dz * dz + dw * dw);
}

