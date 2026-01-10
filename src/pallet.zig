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

    pub fn print(self: *const Self, alloc: std.mem.Allocator) !void {
        var arena = std.heap.ArenaAllocator.init(alloc);
        defer arena.deinit();

        const arena_alloc = arena.allocator();
        _ = &arena_alloc;

        var stdout_writer = std.fs.File.stdout().writer(&.{});
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

