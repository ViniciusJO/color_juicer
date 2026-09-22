const std = @import("std");
const Color = @import("color.zig");

const Self = @This();

primaries: []const *Color.RGBA,
complementaries: []const *Color.RGBA,

pub fn print(self: *const Self, alloc: std.mem.Allocator, io: std.Io) !void {
    var arena = std.heap.ArenaAllocator.init(alloc);
    defer arena.deinit();
    const arena_alloc = arena.allocator();

    var stdout_writer = std.Io.File.stdout().writer(io, &.{});
    const stdout = &stdout_writer.interface;
    var idx: u8 = 1;
    for(self.primaries) |p| {
        defer idx += 1;
        try stdout.print("p{d}: ", .{ idx });
        defer stdout.print("\n", .{}) catch {};
        try print_color(arena_alloc, stdout, p);
    }
    idx = 1;
    try stdout.print("\n", .{});
    for(self.complementaries) |c| {
        defer idx += 1;
        try stdout.print("c{d}: ", .{ idx });
        defer stdout.print("\n", .{}) catch {};
        try print_color(arena_alloc, stdout, c);
    }
}

fn print_color(alloc: std.mem.Allocator, writer: *std.Io.Writer, color: *const Color.RGBA) !void {
    try writer.print("{s} {s}", .{ try color.colorizer(alloc, &"██"),  color.to_rgb_str()  });
}

