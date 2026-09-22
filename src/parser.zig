const std = @import("std");
const Color = @import("color.zig").Color;
const Pallete = @import("pallete.zig");

pub fn capture_to_text(capture: []const u8, pallet: Pallete) []const u8 {
    if(capture.len < 2) return "";
    const number = (std.fmt.parseInt(u8, capture[1..], 10) catch return "") - 1;

    if(capture[0] == 'p') {
        if(number >= pallet.primaries.len) return "";
        return &pallet.primaries[number].to_rgb_str();
    } else if(capture[0] == 'c') {
        if(number >= pallet.complementaries.len) return "";
        return &pallet.complementaries[number].to_rgb_str();
    } else return "";
}

pub fn is_capturable(c: u8) bool {
    return
        (c >= 'a' and c <= 'z') or
        (c >= 'Z' and c <= 'Z') or
        (c >= '0' and c <= '9') or
        (c == '_') or
        (c == '-');
}

pub fn generate_from_template(allocator: std.mem.Allocator, io: std.Io, template: std.Io.File, output: std.Io.File, pallet: Pallete) !void {
    var buff: [4096]u8 = undefined;
    var tr = template.reader(io, &.{});
    var tri = &tr.interface;

    // try tri.readVecAll(&buff);
    const bytes_read = try tri.readSliceShort(&buff);
    // const bytes_read = try template.readAll(&buff);
    const content = buff[0..bytes_read];

    var capture: []const u8 = "";
    var out: []const u8 = "";

    var is_capturing: bool = false;

    for(content, 0..) |c,i| {
        if(is_capturing) {
            if(c == '%') { continue; }
            else if(!is_capturable(c)) {
                const cc = capture_to_text(capture, pallet);
                out = try std.fmt.allocPrint(allocator, "{s}{s}{c}", .{ out, cc, c });
                is_capturing = false;
                capture = "";
                continue;
            }
            capture = try std.fmt.allocPrint(allocator, "{s}{c}", .{ capture, c });
        } else if(c == '%' and ((i < bytes_read - 1) and content[i+1] == c)) {
            is_capturing = true;
        } else out = try std.fmt.allocPrint(allocator, "{s}{c}", .{ out, c });
    }

    var out_reader = output.writer(io, &.{});
    try out_reader.interface.print("{s}", .{ out });
}

pub fn have_sub_path(path: []const u8) bool {
    for(path) |c| {
        if(c == '/') return true;
    }
    return false;
}

pub fn padding(allocator: std.mem.Allocator, level: u8) []const u8 {
    var str: []u8 = "";
    // if(level > 0) str = std.fmt.allocPrint(allocator, "└", .{}) catch str;
    if(level == 0) return str;
    str = std.fmt.allocPrint(allocator, "{s}", .{ str }) catch str;
    for(0..level-1) |_| {
        str = std.fmt.allocPrint(allocator, "{s}\x1b[1;32m│\x1b[0m   ", .{ str }) catch str;
    }

    if(level > 0) str = std.fmt.allocPrint(allocator, "{s}\x1b[1;32m├──\x1b[0m ", .{ str }) catch str
    else str = std.fmt.allocPrint(allocator, "{s}", .{ str }) catch str;
    return str;
}

pub fn iterate_dir_generating_template(
    allocator: std.mem.Allocator,
    io: std.Io,
    ref: std.Io.Dir,
    dir: std.Io.Dir,
    level: u8,
    pallet: Pallete,
) !void {
    var iterable_dir = try dir.walk(allocator);
    while (try iterable_dir.next(io)) |entry| {
        if(have_sub_path(entry.path)) {
            // std.debug.print("<<skip-{s}>>\n", .{entry.path});
            continue;
        }
        std.debug.print("{s}{s} {s}\n", .{
            padding(allocator, level),
            switch(entry.kind) {
                .file => "\x1b[1;33m \x1b[0m",
                .directory => "\x1b[1;34m \x1b[0m",
                else => "  "
            },
            entry.path, 
        });
        // Print the name of each entry
        switch(entry.kind) {
            .file => {
                // std.debug.print("file: {s}\n", .{entry.name});
                const template_file = try dir.openFile(io, entry.path, .{ .mode = .read_only });
                defer template_file.close(io);
                const out_file = ref.openFile(io, entry.path, .{ .mode = .write_only })
                    catch try ref.createFile(io, entry.path, .{});
                defer out_file.close(io);
                try generate_from_template(allocator, io, template_file, out_file, pallet);
            },
            .directory => {
                // std.debug.print("===== DIR =====\n", .{});
                var template_dir = try dir.openDir(io, entry.path, .{ .iterate = true });
                defer template_dir.close(io);
                var out_dir = ref.openDir(io, entry.path, .{ .iterate = true })
                    catch catcher: {
                        try ref.createDir(io, entry.path, .default_dir);
                        // try ref.makeDir(entry.path);
                        break :catcher try ref.openDir(io, entry.path, .{ .iterate = true });
                    };
                defer out_dir.close(io);
                try iterate_dir_generating_template(allocator, io, out_dir, template_dir, level + 1, pallet);
            },
            else => {}
        }
    }
}

pub fn generate_files(io: std.Io, pallet: Pallete) !void {
    const gpa = std.heap.page_allocator;

    var b: [4096]u8 = undefined;

    const home_path = try std.Io.Dir.realPath(io, try std.process.getEnvVarOwned(gpa, "HOME"), &b);
    var home = try std.Io.Dir.openDirAbsolute(io, home_path, .{ .iterate = true });
    defer home.close(io);

    const config_path = try std.fs.path.join(gpa, &[_][]const u8{ home_path, ".config/color_juicer"  });
    std.Io.Dir.createDirAbsolute(io, config_path, .default_dir) catch |e| {
        if(e != error.PathAlreadyExists) return e;
    };

    const template_path = try std.fs.path.join(gpa, &[_][]const u8{ config_path, "template"  });
    std.Io.Dir.createDirAbsolute(io, template_path, .default_dir) catch |e| {
        if(e != error.PathAlreadyExists) return e;
    };
    var template = try std.Io.Dir.openDirAbsolute(io, template_path, .{ .iterate = true });
    defer template.close(io);

    try iterate_dir_generating_template(gpa, io, home, template, 0, pallet);
}

