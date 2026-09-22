const std = @import("std");
const args_parser = @import("args_parser");

const stb_image = @import("stbi");
const stb_image_write = @import("stbiw");
const stb_image_resize = @import("stbir");

const Color = @import("color.zig");
const Vec = @import("vec.zig");
const Pallet = @import("pallet.zig").Pallet;
const Pallet_ = @import("pallet.zig").Pallet_;
const Parser = @import("parser.zig");

const write_downsampled_image = false;

const Image = struct {
    width: c_int,
    height: c_int,
    channels: c_int,
    pixels: []u8,
};

const Mean = struct {
    const Self = @This();

    color: Vec.Vec4,
    colors: std.ArrayList(Vec.Vec4),
    dist: f64,
    partition_size: u64 = 0,

    pub fn print_color(self: *const Self, alloc: std.mem.Allocator, io: std.Io) !void {
        var arena = std.heap.ArenaAllocator.init(alloc);
        defer arena.deinit();
        const color = Color.RGBA.init(.{ .vec = self.color });
        const arena_alloc = arena.allocator();
        _ = &arena_alloc;
        var stdout_writer = std.Io.File.stdout().writer(io, &.{});
        const stdout = &stdout_writer.interface;
        try stdout.print("{s} {s}\n", .{ try color.colorizer(arena_alloc, &"██"), color.to_rgb_str() });
        try stdout.flush();
    }
};

pub fn kmeanspp_init(alloc: std.mem.Allocator, io: std.Io, m: *[]Mean, pixels: *[]const Vec.Vec3) !void {
    const n = pixels.len;
    const k = m.len;

    const first_index = try randomIndex(io, n);
    const first_pixel = pixels.*[first_index];

    m.*[0].color = .{ first_pixel[0], first_pixel[1], first_pixel[2], 255 };

    var ar = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer ar.deinit();
    std.debug.print("random_i: {} {any} {any} {s}\n", .{
        first_index,
        first_pixel,
        m.*[0].color,
        Color.RGBA.init(.{ .vec = m.*[0].color }).to_rgb_str()
    });

    var distances = try alloc.alloc(f64, n);
    defer alloc.free(distances);

    for (1..k) |i| {
        for (0..pixels.len) |j| {
            // std.debug.print("pixel: {}\n", .{j});
            // 221178
            const p = pixels.*[j];
            var min_dist: f64 = std.math.inf(f64);
            for (m.*[0..i]) |mean| {
                const dist_sq = vec_dist_sq(mean.color, .{ p[0], p[1], p[2], 255 });
                if (dist_sq < min_dist) min_dist = dist_sq;
            }
            distances[j] = min_dist;
        }

        var total: f64 = 0;
        for (distances) |d| total += d;

        const r = try randomFloat(io) * total;
        var cumulative: f64 = 0;
        var next_index: usize = 0;
        for (distances, 0..) |d, j| {
            cumulative += d;
            if (cumulative >= r) {
                next_index = j;
                break;
            }
        }
        const np = pixels.*[next_index];
        m.*[i].color = .{ np[0], np[1], np[2], 255 };
    }
}

fn repartition(alloc: std.mem.Allocator, ms: *[]Mean, pixels: *[]const Vec.Vec3) !*[]Mean {
    for (ms.*) |*m| m.*.colors.clearRetainingCapacity();
    for (pixels.*) |p| {
        // pixels are Vec3 (matches the image's 3 channels), Mean.color is
        // Vec4 -- pad with a fixed alpha, same convention as kmeanspp_init.
        const c: Vec.Vec4 = .{ p[0], p[1], p[2], 255 };
        var min = &(ms.*[0]);
        for (ms.*) |*m| {
            if (Vec.vec4_dist(c, m.color) < Vec.vec4_dist(c, min.color)) min = m;
        }
        try min.colors.append(alloc, c);
        min.partition_size +|= 1;
    }
    return ms;
}

fn compute_means(ms: *[]Mean, min_dist: f64) bool {
    var updated: usize = 0;
    for (ms.*) |*m| {
        if (m.*.dist <= min_dist) continue;

        var sum: [4]usize = .{0} ** 4;
        for (m.*.colors.items) |color| {
            sum[0] += color[0];
            sum[1] += color[1];
            sum[2] += color[2];
            sum[3] += color[3];
        }

        const color_res = if(m.*.colors.items.len > 0) Vec.Vec4 {
            @intCast(sum[0] / m.*.colors.items.len),
            @intCast(sum[1] / m.*.colors.items.len),
            @intCast(sum[2] / m.*.colors.items.len),
            @intCast(sum[3] / m.*.colors.items.len)
        } else Vec.Vec4{ 0 , 0, 0, 0 };

        m.*.dist = Vec.vec4_dist(m.*.color, color_res);
        m.*.color = color_res;
        updated += 1;
    }
    return updated == 0;
}

fn vec_dist_sq(a: Vec.Vec4, b: Vec.Vec4) f64 {
    const square = struct { pub fn sq(x: f64) f64 { return x*x; } }.sq;
    return
        square(@as(f64, @floatFromInt(a[0])) - @as(f64, @floatFromInt(b[0]))) + 
        square(@as(f64, @floatFromInt(a[1])) - @as(f64, @floatFromInt(b[1]))) + 
        square(@as(f64, @floatFromInt(a[2])) - @as(f64, @floatFromInt(b[2]))) + 
        square(@as(f64, @floatFromInt(a[3])) - @as(f64, @floatFromInt(b[3])));
}

fn randomIndex(io: std.Io, n: usize) !usize {
    var buf: u64 = undefined;
    // std.crypto.random.bytes(std.mem.asBytes(&buf));
    std.Io.random(io, std.mem.asBytes(&buf));
    return @intCast(buf % @as(u64, n));
}

fn randomFloat(io: std.Io) !f64 {
    var buf: u64 = undefined;
    // std.crypto.random.bytes(std.mem.asBytes(&buf));
    std.Io.random(io, std.mem.asBytes(&buf));
    // Divide by max u64 to get float in [0,1)
    return @as(f64, @floatFromInt(buf)) / @as(f64, @floatFromInt(std.math.maxInt(u64)));
}

fn usage(name: []const u8) void {
    // std.debug.print("Usage {s}:\n\t{s} filepath <#means> <downsampling_factor>\n", .{name, name});
    std.debug.print("Usage:\t{s} filepath\n", .{name});
}

const Args = struct {
    image_path: []const u8,
    contrast: ?[]const u8,
    templates_path: ?[]const u8,
    output_reference_path: ?[]const u8,
};

const StdIO = struct {
    in: *std.Io.Reader,
    out: *std.Io.Writer,
    err: *std.Io.Writer,
};

const Kit = struct {
    const Self = @This();

    io: std.Io,
    args: Args,
    stdio: StdIO,
    gpa: std.mem.Allocator,
    arena: std.heap.ArenaAllocator,

    pub fn init(io: std.Io, arguments: Args) Self {
        var stdin_reader = std.Io.File.stdin().reader(io, &.{});
        var stdout_writer = std.Io.File.stdout().writer(io, &.{});
        var stderr_writer = std.Io.File.stderr().writer(io, &.{});

        const stdin = &stdin_reader.interface;
        const stdout = &stdout_writer.interface;
        const stderr = &stderr_writer.interface;

        return .{
            .io = io,
            .args = arguments,
            .gpa = std.heap.page_allocator,
            .arena = std.heap.ArenaAllocator.init(std.heap.page_allocator),
            .stdio = StdIO { .in = stdin, .out = stdout, .err = stderr, }
        };
    }
};

pub fn main(init: std.process.Init) !void {

    const alloc = init.gpa;
    const io = init.io;

    // var stdin_reader = std.Io.File.stdin().reader(io, &.{});
    // const stdin = &stdin_reader.interface;

    var stdout_writer = std.Io.File.stdout().writer(io, &.{});
    const stdout = &stdout_writer.interface;

    var stderr_writer = std.Io.File.stderr().writer(io, &.{});
    const stderr = &stderr_writer.interface;

    const args = try init.minimal.args.toSlice(alloc);
    defer alloc.free(args);
    errdefer usage(args[0][0..args[0].len]);

    var _a: Args = undefined;
    _a.image_path = if(args.len >= 1) args[1] else {
        try stderr.print("ERROR: missing required file path\n", .{});
        try stderr.flush();
        return;
    };

    // const means_quant = if(args.len >= 2) std.fmt.parseInt(usize, args[2], 10) catch return error.InvalidArgument else 10;
    // const fac = if(args.len >= 3) std.fmt.parseInt(u8, args[3], 10) catch return error.InvalidArgument else 4;

    const arguments = _a;

    // const kit = Kit{
    //     .io = io,
    //     .args = arguments,
    //     .gpa = alloc,
    //     .arena = std.heap.ArenaAllocator.init(alloc),
    //     .stdio = StdIO { .in = stdin, .out = stdout, .err = stderr, }
    // };

    // const kit = Kit.init(io,arguments);
    // _ = kit;


    // const input = try kit.stdio.in.takeDelimiterExclusive('\n');//readSliceShort(&bff);
    //
    // try kit.stdio.out.print("\nTEST: {s}\n", .{ input });

    // if(true) return;




    const means_quant = 6;
    const fac = 5;
    const precision = 0.01;
    // const precision = 1;
    // const precision = 100;
    // const precision = 1000;
    // #112313
    
    var input_image: Image = undefined;

    const input_c_ptr = stb_image.stbi_load(
        arguments.image_path.ptr,
        &input_image.width,
        &input_image.height,
        &input_image.channels,
        0
    );
    if(null == input_c_ptr) return error.ImageLoadError;
    input_image.pixels = input_c_ptr[0..@as(usize, @intCast(input_image.width*input_image.height*input_image.channels))];
    defer std.c.free(input_c_ptr);
    var downsampled = Image{
        .width = @divTrunc(input_image.width, fac),
        .height = @divTrunc(input_image.height, fac),
        .channels = 3,
        .pixels = undefined,
    };
    
    const down_c_ptr = stb_image_resize.stbir_resize_uint8_linear(
        input_image.pixels.ptr,
        input_image.width, input_image.height, input_image.width*input_image.channels*@sizeOf(u8),
        null,
        downsampled.width, downsampled.height, downsampled.width*downsampled.channels*@sizeOf(u8),
        @intCast(downsampled.channels)
    );
    if(null == down_c_ptr) return error.ImageDownsamplingError;
    defer std.c.free(down_c_ptr.?);
    downsampled.pixels = down_c_ptr[0..@as(usize, @intCast(downsampled.width*downsampled.height*downsampled.channels))];

    std.debug.print("\nInput: Image{{ .width = {}, .height = {}, .channels = {} }}\n Down: Image{{ .width = {}, .height = {}, .channels = {} }}\n\n", .{ input_image.width, input_image.height, input_image.channels, downsampled.width, downsampled.height, downsampled.channels });
    std.debug.print("down_length: {} bytes\n", .{ downsampled.pixels.len });

    if(write_downsampled_image) _ = stb_image_write.stbi_write_png(
        "dout.png",
        downsampled.width,
        downsampled.height,
        downsampled.channels,
        downsampled.pixels.ptr,
        downsampled.width*downsampled.channels*@sizeOf(u8)
    );

    var means = try alloc.alloc(Mean, means_quant);
    defer alloc.free(means);

    for(means) |*m| {
        m.color = Color.RGBA.init(.{ .hex = 0x000000FF }).to_vec4();
        m.colors = try std.ArrayList(Vec.Vec4).initCapacity(alloc, 0);
        m.dist = std.math.inf(f64);
        m.partition_size = 0;
    }

    var arena = std.heap.ArenaAllocator.init(alloc);
    defer arena.deinit();

    var pixel_view: []const Vec.Vec3 = Vec.bytes_as_vec3(downsampled.pixels);

    try kmeanspp_init(arena.allocator(), io, &means, &pixel_view);

    while (true) {
        _ = try repartition(arena.allocator(), &means, &pixel_view);
        if (compute_means(&means, precision)) break;
    }

    for (means) |*m| {
        if(m.*.partition_size == 0) m.*.color = means[0].color;
    }

    // std.debug.print("\nMEANS:\n", .{});
    for (means) |m| try m.print_color(alloc, io);
    // std.debug.print("\n\n", .{});

    const color_dist_square = struct { fn color_dist(_clr: Color.RGBA) u64 {
        const clr = [_]u64{ _clr.r, _clr.g, _clr.b, _clr.a };
        return clr[0]*clr[0] + clr[1]*clr[1] + clr[2]*clr[2];
    }}.color_dist;

    std.sort.heap(Mean, means, {}, struct {
        pub fn cmp(_: void, a: Mean, b: Mean) bool {
            const __a = Color.RGBA.init(.{ .vec = a.color });
            const __b = Color.RGBA.init(.{ .vec = b.color });
            const c1_ = color_dist_square(__a);
            const c2_ = color_dist_square(__b);

            const a_ = __a.to_lch();
            const b_ = __b.to_lch();
            const c1 = a_.c/a_.l*@as(f32, @floatFromInt(c1_));
            const c2 = b_.c/b_.l*@as(f32, @floatFromInt(c2_));

            return
                a.partition_size > b.partition_size or
                c1 > c2;
        }
    }.cmp);
    
    var primary = try alloc.alloc(Color.RGBA, means.len);
    defer alloc.free(primary);
    var complementary = try alloc.alloc(Color.RGBA, means.len);
    defer alloc.free(complementary);

    for (means, 0..) |*m, i| {
        primary[i] = Color.RGBA.init(.{ .vec = m.color });
        complementary[i] = primary[i].complementary();
    }

    //
    // for (primary) |*c| {
    //     if(c.to_lch().l < 0.5) {
    //         var aux = c.to_lch();
    //         aux.l = 0.6;
    //         c.* = aux.to_rgba();
    //     }
    // }

    const pallet = Pallet{
        .prim  = &primary[0],
        .sec   = &primary[1],
        .terc  = &primary[2],
        .cprim = &complementary[0],
        .csec  = &complementary[1],
        .cterc = &complementary[2],
    };

    const pallet_ = Pallet_{
        .primaries = &[_]*Color.RGBA{ &primary[0], &primary[1], &primary[2], },
        .complementaries = &[_]*Color.RGBA{ &complementary[0],  &complementary[1], &complementary[2], },
    };

    try stdout.print("\nPallet:\n\n", .{});
    try pallet.print(alloc, io);

    try stdout.print("\nPallet_:\n\n", .{});
    try pallet_.print(alloc, io);


    // TODO: get input directory from args
    const in_path  = "src/.ignore/template";
    const out_path = "src/.ignore/out";

    var realpath_in_buf: [4096]u8 = undefined;
    const inx = std.Io.Dir.realPathFile(std.Io.Dir.cwd(), init.io, in_path, &realpath_in_buf) catch |err| {
        std.log.err("Failed to open path \"{s}\": {}", .{ in_path, err });
        std.process.exit(1);
    };
    // TODO: verify if the realpath points to actual dirs
    try stdout.print("\nIN__REALPAT: {s}\n", .{realpath_in_buf[0..inx]});
    var in  = try std.Io.Dir.openDirAbsolute(io, realpath_in_buf[0..inx], .{ .iterate = true, .access_sub_paths = true });
    defer in.close(io);

    var realpath_out_buf: [4096]u8 = undefined;
    const outx = std.Io.Dir.realPathFile(std.Io.Dir.cwd(), init.io, out_path, &realpath_out_buf) catch |err| {
        std.log.err("Failed to open path \"{s}\": {}", .{ out_path, err });
        std.process.exit(1);
    };
    // TODO: verify if the realpath points to actual dirs
    try stdout.print("OUT_REALPAT: {s}\n\n", .{realpath_out_buf[0..outx]});
    var out = try std.Io.Dir.openDirAbsolute(io, realpath_out_buf[0..outx], .{ .iterate = true, .access_sub_paths = true });
    defer out.close(io);


    try stdout.print("Generating: files from templates...\n\n", .{});
    try Parser.iterate_dir_generating_template(arena.allocator(), io, out, in, 0, pallet);






    // var b: [4096]u8 = undefined;
    //
    // const home_path = try std.fs.realpath(try std.process.getEnvVarOwned(gpa, "HOME"), &b);
    // var home = try std.fs.openDirAbsolute(home_path, .{ .iterate = true });
    // defer home.close();
    //
    // const config_path = try std.fs.path.join(gpa, &[_][]const u8{ home_path, ".config/color_juicer"  });
    // std.fs.makeDirAbsolute(config_path) catch |e| {
    //     if(e != error.PathAlreadyExists) return e;
    // };
    //
    // const template_path = try std.fs.path.join(gpa, &[_][]const u8{ config_path, "template"  });
    // std.fs.makeDirAbsolute(template_path) catch |e| {
    //     if(e != error.PathAlreadyExists) return e;
    // };
    // var template = try std.fs.openDirAbsolute(template_path, .{ .iterate = true });
    // defer template.close();







    // ========================= ARCHIVE =========================


    // for (means) |*m| {
    //     std.debug.print("{any}\n", .{m.*.color});
    //     std.debug.print("#{x}{x}{x}{x}\n", .{m.*.color[0], m.*.color[1], m.*.color[2], m.*.color[3]});
    //     // std.debug.print("\n{s}\n", .{try color_string(&"██", m.*.color)});
    //     try stdout.print("{} {} {}\n", .{m.*.color[0], m.*.color[1], m.*.color[2]});
    // }





    // var comp_mean: [3]u64 = .{0, 0, 0};
    // var sum: u64 = 1;
    // for (complementary) |*m| {
    //     // sum += m.*.partition_size;
    //     comp_mean[0] += m.r * m.*.partition_size;
    //     comp_mean[1] += m.g * m.*.partition_size;
    //     comp_mean[2] += m.b * m.*.partition_size;
    // }
    // comp_mean[0] /= sum;
    // comp_mean[1] /= sum;
    // comp_mean[2] /= sum;

    
    // const clrs = means;
    // const cclrs = complementary;


    // try std.io.getStdOut().writer().print("prim: {s}\n", .{ try color_string_(&"██", clrs[0].color) });
    // try std.io.getStdOut().writer().print("sec: {s}\n", .{ try color_string_(&"██", clrs[1].color) });
    // try std.io.getStdOut().writer().print("cprim: {s}\n", .{ try color_string_(&"██", cclrs[0].color) });
    // try std.io.getStdOut().writer().print("csec: {s}\n", .{ try color_string_(&"██", cclrs[1].color) });
    // try std.io.getStdOut().writer().print("cont: {s}\n", .{ try color_string_(&"██", comp_m) });

    // const color_dist_square = struct { fn color_dist(_clr: Color.RGBA) u64 {
    //     const clr = [_]u64{ _clr.r, _clr.g, _clr.b, _clr.a };
    //     return clr[0]*clr[0] + clr[1]*clr[1] + clr[2]*clr[2];
    // }}.color_dist;

    // std.sort.heap(Mean, means, {}, struct {
    //     pub fn cmp(_: void, _a: Mean, _b: Mean) bool {
    //         const a = _a.color;
    //         const b = _b.color;
    //         const c1_ = color_dist_square(_a.color);
    //         const c2_ = color_dist_square(_b.color);
    //
    //         const a_ = a.to_lch();
    //         const b_ = b.to_lch();
    //         const c1 = a_.c*a_.l*@as(f32, @floatFromInt(c1_));
    //         const c2 = b_.c*b_.l*@as(f32, @floatFromInt(c2_));
    //
    //         return c1 > c2;
    //     }
    // }.cmp);

    // for(0..4) |_| {
    //     for (means.items) |*m| {
    //         const col = m.*.color;
    //         std.debug.print("{s}", .{try color_string(&"████████████", col)});
    //     }
    //     std.debug.print("\n", .{});
    // }
    // for (means.items) |*m| {
    //     std.debug.print("{d:^12}", .{ m.*.partition_size });
    // }
    // std.debug.print("\n\n", .{});

    // {
    //     // var out_writer = std.fs.File.stdout().writer(&buff).interface;
    //     var o = out_ini_file.writer(&.{});
    //     var out_out_1 = &o.interface;
    //     try out_out_1.print("pprim = #{X}{X}{X}\n", .{ clrs[0].color[0],  clrs[0].color[1], clrs[0].color[2] });
    //     try out_out_1.print("psec = #{X}{X}{X}\n", .{ clrs[1].color[0],  clrs[1].color[1], clrs[1].color[2] });
    //     try out_out_1.print("pterc = #{X}{X}{X}\n", .{ clrs[2].color[0],  clrs[2].color[1], clrs[2].color[2] });
    //     try out_out_1.print("pcont = #{X}{X}{X}\n", .{ comp_m[0], comp_m[1], comp_m[2] });
    //     // try out_writer.flush();
    // }

    // try out_i3_file.writer().print("set $text_focus   #{X}{X}{X}\n", .{ clrs[0].color[0],  clrs[0].color[1], clrs[0].color[2] });
    // try out_i3_file.writer().print("set $bg_normal    #{X}{X}{X}\n", .{ clrs[1].color[0],  clrs[1].color[1], clrs[1].color[2] });
    // try out_i3_file.writer().print("set $text_normal  #{X}{X}{X}\n", .{ clrs[2].color[0],  clrs[2].color[1], clrs[2].color[2] });
    // try out_i3_file.writer().print("set $bg_focus     #{X}{X}{X}\n", .{ comp_m[0], comp_m[1], comp_m[2] });
    // borda | fundo título | texto título | indicador | texto título (estado inverso)
    // try out_i3_file.writer().print("\nclient.focused    $bg_focus $bg_focus #000000 $bg_focus $bg_focus\n", .{});

    // const color_to_rgb_str = struct {
    //     fn col(co: Color.RGBA) [7]u8 {
    //         var ret: [7]u8 = undefined;
    //         var stream = std.io.fixedBufferStream(&ret);
    //         stream.writer().print("#{X:02}{X:02}{X:02}", .{ co.r, co.g, co.b }) catch { ret = .{ '#', '0', '0', '0', '0', '0', '0' }; };
    //         return ret;
    //     }
    // }.col;

    // const Colorscheme = struct {
    //     primary: [7]u8,
    //     secondary: [7]u8,
    //     terciary: [7]u8,
    //     complementary: [7]u8
    // };

    // const cs = Colorscheme {
    //     .primary  = color_to_rgb_str(clrs[0].color),
    //     .secondary  = color_to_rgb_str(clrs[0].color),
    //     .terciary  = color_to_rgb_str(clrs[0].color),
    //     .complementary  = color_to_rgb_str(Color.RGBA.init(.{ .vec = comp_m })),
    // };


    // TODO: I3
    // // borda | fundo título | texto título | indicador | texto título (estado inverso)
    // try out_i3_writeri.print("\nclient.focused {s} {s} {s} {s} {s}\n", .{
    //     cs.complementary, // title border
    //     cs.complementary, // title background
    //     "#000000",        // title text
    //     cs.terciary     , // indicator
    //     cs.complementary  // border
    // });




    // TODO: polybar
    //
    // // [dyn_colors]
    // // pprim = #D8ACBB
    // // psec = #B68DAA
    // // pterc = #9C769E
    // // pcont = #5C6C48
    // // cprim = #496048
    // // csec = #627B55
    // 
        // try out_writeri.print("[dyn_colors]\n", .{});
        // try out_writeri.print("pprim = #{X}{X}{X}\n", .{ clrs[0].color.r,  clrs[0].color.g, clrs[0].color.b });
        // try out_writeri.print("psec = #{X}{X}{X}\n", .{ clrs[1].color.r,  clrs[1].color.g, clrs[1].color.b });
        // try out_writeri.print("pterc = #{X}{X}{X}\n", .{ clrs[2].color.r,  clrs[2].color.g, clrs[2].color.b });
        // try out_writeri.print("pcont = #{X}{X}{X}\n", .{ comp_m[0], comp_m[1], comp_m[2] });
        //
        // try out_writeri.print("cprim = #{X}{X}{X}\n", .{ cclrs[0].color.r,  cclrs[0].color.g, cclrs[0].color.b });
        // try out_writeri.print("csec = #{X}{X}{X}\n", .{ cclrs[1].color.r,  cclrs[1].color.g, cclrs[1].color.b });
        // try out_writeri.print("cont = #{X}{X}{X}\n", .{ comp_m[0], comp_m[1], comp_m[2] });

        // try out_out_1.print("pprim = #{X}{X}{X}\n", .{ clrs[0].color[0],  clrs[0].color[1], clrs[0].color[2] });
        // try out_out_1.print("psec = #{X}{X}{X}\n", .{ clrs[1].color[0],  clrs[1].color[1], clrs[1].color[2] });
        // try out_out_1.print("pterc = #{X}{X}{X}\n", .{ clrs[2].color[0],  clrs[2].color[1], clrs[2].color[2] });
        // try out_out_1.print("pcont = #{X}{X}{X}\n", .{ comp_m[0], comp_m[1], comp_m[2] });

        // std.debug.print("[dyn_colors]\n", .{});
        // std.debug.print("prim = #{X}{X}{X}\n", .{ clrs[0].color[0],  clrs[0].color[1], clrs[0].color[2] });
        // std.debug.print("sec = #{X}{X}{X}\n", .{ clrs[1].color[0],  clrs[1].color[1], clrs[1].color[2] });
        // std.debug.print("cprim = #{X}{X}{X}\n", .{ cclrs[0].color[0],  cclrs[0].color[1], cclrs[0].color[2] });
        // std.debug.print("csec = #{X}{X}{X}\n", .{ cclrs[1].color[0],  cclrs[1].color[1], cclrs[1].color[2] });
        // std.debug.print("cont = #{X}{X}{X}\n", .{ comp_m[0], comp_m[1], comp_m[2] });
        // try out_writer.flush();
    // }
}
