const std = @import("std");
const Color = @import("color.zig").Color;
const Pallet = @import("pallet.zig").Pallet;

pub fn capture_to_text(capture: []const u8, pallet: Pallet) []const u8 {

    // if(std.mem.eql(u8, "prim", capture))      { return &pallet.prim.rgb_str(); }
    // else if(std.mem.eql(u8, "sec", capture))  { return &pallet.sec.rgb_str(); }
    // else if(std.mem.eql(u8, "terc", capture)) { return &pallet.terc.rgb_str(); }
    // else if(std.mem.eql(u8, "comp", capture)) { return &pallet.comp.rgb_str(); }
    // else return "";
    
    if(std.mem.eql(u8, "prim", capture))       { return &pallet.prim.to_rgb_str();  }
    else if(std.mem.eql(u8, "sec", capture))   { return &pallet.sec.to_rgb_str();   }
    else if(std.mem.eql(u8, "terc", capture))  { return &pallet.terc.to_rgb_str();  }
    else if(std.mem.eql(u8, "cprim", capture)) { return &pallet.cprim.to_rgb_str(); }
    else if(std.mem.eql(u8, "csec", capture))  { return &pallet.csec.to_rgb_str();  }
    else if(std.mem.eql(u8, "cterc", capture)) { return &pallet.cterc.to_rgb_str(); }
    else return "";
}

pub fn is_capturable(c: u8) bool {
    return
        (c >= 'a' and c <= 'z') or
        (c >= 'Z' and c <= 'Z') or
        (c >= '0' and c <= '9') or
        (c == '_') or
        (c == '-');
}

pub fn generate_from_template(allocator: std.mem.Allocator, template: std.fs.File, output: std.fs.File, pallet: Pallet) !void {
    var buff: [4096]u8 = undefined;
    const bytes_read = try template.readAll(&buff);
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

    var out_reader = output.writer(&.{});
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
    ref: std.fs.Dir,
    dir: std.fs.Dir,
    level: u8,
    pallet: Pallet,
) !void {
    var iterable_dir = try dir.walk(allocator);
    while (try iterable_dir.next()) |entry| {
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
                const template_file = try dir.openFile(entry.path, .{ .mode = .read_only });
                defer template_file.close();
                const out_file = ref.openFile(entry.path, .{ .mode = .write_only })
                    catch try ref.createFile(entry.path, .{});
                defer out_file.close();
                try generate_from_template(allocator, template_file, out_file, pallet);
            },
            .directory => {
                // std.debug.print("===== DIR =====\n", .{});
                var template_dir = try dir.openDir(entry.path, .{ .iterate = true });
                defer template_dir.close();
                var out_dir = ref.openDir(entry.path, .{ .iterate = true })
                    catch catcher: {
                        try ref.makeDir(entry.path);
                        break :catcher try ref.openDir(entry.path, .{ .iterate = true });
                    };
                defer out_dir.close();
                try iterate_dir_generating_template(allocator, out_dir, template_dir, level + 1, pallet);
            },
            else => {}
        }
    }
}

pub fn generate_files(pallet: Pallet) !void {
    const gpa = std.heap.page_allocator;

    var b: [4096]u8 = undefined;

    const home_path = try std.fs.realpath(try std.process.getEnvVarOwned(gpa, "HOME"), &b);
    var home = try std.fs.openDirAbsolute(home_path, .{ .iterate = true });
    defer home.close();

    const config_path = try std.fs.path.join(gpa, &[_][]const u8{ home_path, ".config/color_juicer"  });
    std.fs.makeDirAbsolute(config_path) catch |e| {
        if(e != error.PathAlreadyExists) return e;
    };

    const template_path = try std.fs.path.join(gpa, &[_][]const u8{ config_path, "template"  });
    std.fs.makeDirAbsolute(template_path) catch |e| {
        if(e != error.PathAlreadyExists) return e;
    };
    var template = try std.fs.openDirAbsolute(template_path, .{ .iterate = true });
    defer template.close();

    try iterate_dir_generating_template(gpa, home, template, 0, pallet);
}

