const std = @import("std");
const Color = @import("color.zig");

pub const Pallet = struct {
    prim:  *Color.RGBA,
    sec:   *Color.RGBA,
    terc:  *Color.RGBA,
    cprim: *Color.RGBA,
    csec:  *Color.RGBA,
    cterc: *Color.RGBA,

    const Self = @This();

    pub fn print(self: *const Self, alloc: std.mem.Allocator, io: std.Io) !void {
        var arena = std.heap.ArenaAllocator.init(alloc);
        defer arena.deinit();

        const arena_alloc = arena.allocator();
        _ = &arena_alloc;


        // var stdout_writer = std.fs.File.stdout().writer(&.{});
        var stdout_writer = std.Io.File.stdout().writer(io, &.{});
        const stdout = &stdout_writer.interface;
        try stdout.print("prim:  {s} {s}\n", .{ try self.prim.colorizer(arena_alloc, &"██"),  self.prim.to_rgb_str()  });
        try stdout.print("sec:   {s} {s}\n", .{ try self.sec.colorizer(arena_alloc, &"██"),   self.sec.to_rgb_str()   });
        try stdout.print("terc:  {s} {s}\n", .{ try self.terc.colorizer(arena_alloc, &"██"),  self.terc.to_rgb_str()  });
        try stdout.print("cprim: {s} {s}\n", .{ try self.cprim.colorizer(arena_alloc, &"██"), self.cprim.to_rgb_str() });
        try stdout.print("csec:  {s} {s}\n", .{ try self.csec.colorizer(arena_alloc, &"██"),  self.csec.to_rgb_str()  });
        try stdout.print("cterc: {s} {s}\n", .{ try self.cterc.colorizer(arena_alloc, &"██"), self.cterc.to_rgb_str() });
        try stdout.flush();
    }
};

fn print_color(alloc: std.mem.Allocator, writer: *std.Io.Writer, color: *const Color.RGBA) !void {
    try writer.print("{s} {s}", .{ try color.colorizer(alloc, &"██"),  color.to_rgb_str()  });
}

pub const Pallet_ = struct {
    primaries: []const *Color.RGBA,
    complementaries: []const *Color.RGBA,

    const Self = @This();

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
};

