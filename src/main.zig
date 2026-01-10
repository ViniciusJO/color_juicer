const std = @import("std");
const Color = @import("color.zig");
const Vec = @import("vec.zig");
const Pallet = @import("pallet.zig").Pallet;
const Parser = @import("parser.zig");

const stb_image = @import("stbi");
// const stb_image_write = @import("stbiw");
const stb_image_resize = @import("stbir");

const Mean = struct { color: Vec.Vec4, colors: std.ArrayList(Vec.Vec4), dist: f64, partition_size: u64 = 0 };

pub fn kmeanspp_init(alloc: std.mem.Allocator, m: *[]Mean, pixels: *[]const Vec.Vec4) !void {
    const n = pixels.len;
    const k = m.len;

    // Pick first centroid randomly
    const first_index = try randomIndex(n);
    const first_pixel = pixels.*[first_index];

    m.*[0].color = first_pixel;

    // Allocate distances array
    var distances = try alloc.alloc(f64, n);
    defer alloc.free(distances);

    // For each remaining centroid
    for (1..k) |i| {
        // For each pixel, find distance squared to nearest existing centroid
        for (pixels.*, 0..) |p, j| {
            var min_dist: f64 = std.math.inf(f64);
            for (m.*[0..i]) |mean| {
                const dist_sq = vec_dist_sq(mean.color, p);
                if (dist_sq < min_dist)
                    min_dist = dist_sq;
            }
            distances[j] = min_dist;
        }

        // Calculate total weighted distance
        var total: f64 = 0;
        for (distances) |d| total += d;

        // Select next centroid with probability proportional to distance squared
        const r = try randomFloat() * total;
        var cumulative: f64 = 0;
        var next_index: usize = 0;
        for (distances, 0..) |d, j| {
            cumulative += d;
            if (cumulative >= r) {
                next_index = j;
                break;
            }
        }

        // Set the new centroid
        m.*[i].color = pixels.*[next_index];
    }
}

fn repartition(alloc: std.mem.Allocator, ms: *[]Mean, pixels: *[]const Vec.Vec4) !*[]Mean {
    for (ms.*) |*m| m.*.colors.clearRetainingCapacity(); //.clearRetainingCapacity();
    for (pixels.*) |c| {
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

fn randomIndex(n: usize) !usize {
    var buf: u64 = undefined;
    std.crypto.random.bytes(std.mem.asBytes(&buf));
    return @intCast(buf % @as(u64, n));
}

fn randomFloat() !f64 {
    var buf: u64 = undefined;
    std.crypto.random.bytes(std.mem.asBytes(&buf));
    // Divide by max u64 to get float in [0,1)
    return @as(f64, @floatFromInt(buf)) / @as(f64, @floatFromInt(std.math.maxInt(u64)));
}

fn usage(name: []u8) void {
    // std.debug.print("Usage {s}:\n\t{s} filepath <#means> <downsampling_factor>\n", .{name, name});
    std.debug.print("Usage:\t{s} filepath\n", .{name});
}

const Args = struct {
    image_path: []const u8,
    templates_path: []const u8,
    output_reference_path: []const u8,
};

pub fn main() !void {
    const gpa = std.heap.page_allocator;

    const args = try std.process.argsAlloc(gpa);
    defer std.process.argsFree(gpa, args);
    errdefer usage(args[0]);

    var _a: Args = undefined;
    _a.image_path = if(args.len >= 1) args[1] else return error.NoFilepath;
    // const means_quant = if(args.len >= 2) std.fmt.parseInt(usize, args[2], 10) catch return error.InvalidArgument else 10;
    // const fac = if(args.len >= 3) std.fmt.parseInt(u8, args[3], 10) catch return error.InvalidArgument else 4;
    const arguments = _a;

    const means_quant = 3;
    const fac = 5;
    const precision = 0.01;
    

    const Image = struct {
        width: c_int,
        height: c_int,
        channels: c_int,
        pixels: []u8,
    };

    var input_image: Image = undefined;

    input_image.pixels = (stb_image.stbi_load(
        arguments.image_path.ptr,
        &input_image.width,
        &input_image.height,
        &input_image.channels,
        0
    ))[0..@as(usize, @intCast(input_image.width*input_image.height*input_image.channels))];
    // if(null == input_image.pixels) return error.ImageLoadError;

    var downsampled = Image{
        .width = @divTrunc(input_image.width, fac),
        .height = @divTrunc(input_image.height, fac),
        .channels = 4,
        .pixels = undefined,
    };
    
    downsampled.pixels = (stb_image_resize.stbir_resize_uint8_linear(
        input_image.pixels.ptr,
        input_image.width, input_image.height, input_image.width*input_image.channels*@sizeOf(u8),
        null,
        downsampled.width, downsampled.height, downsampled.width*downsampled.channels*@sizeOf(u8),
        @intCast(downsampled.channels)
    ))[0..@as(usize, @intCast(downsampled.width*downsampled.height*downsampled.channels))];
    // if(null == downsampled) return error.ImageDownsamplingError;

    // const down_size: usize = @intCast(new_w*new_h*new_c);
    // var down_slice: []const u8 = downsampled[0..down_size];

    var means = try gpa.alloc(Mean, means_quant);
    defer gpa.free(means);

    for(means) |*m| {
        m.color = Color.RGBA.init(.{ .hex = 0x000000FF }).to_vec4();
        m.colors = try std.ArrayList(Vec.Vec4).initCapacity(gpa, 0);
        m.dist = std.math.inf(f64);
        m.partition_size = 0;
    }
    try kmeanspp_init(gpa, &means, @ptrCast(&downsampled.pixels));

    while (true) {
        _ = try repartition(gpa, &means, @ptrCast(&downsampled.pixels));
        if (compute_means(&means, precision)) break;
    }

    for (means) |*m| {
        if(m.*.partition_size == 0) m.*.color = means[0].color;
    }

    // const color_dist_square = struct { fn color_dist(_clr: Color.RGBA) u64 {
    //     const clr = [_]u64{ _clr.r, _clr.g, _clr.b, _clr.a };
    //     return clr[0]*clr[0] + clr[1]*clr[1] + clr[2]*clr[2];
    // }}.color_dist;

    // std.sort.heap(Mean, means, {}, struct {
    //     pub fn cmp(_: void, a: Mean, b: Mean) bool {
    //         const __a = Color.RGBA.init(.{ .vec = a.color });
    //         const __b = Color.RGBA.init(.{ .vec = b.color });
    //         const c1_ = color_dist_square(__a);
    //         const c2_ = color_dist_square(__b);
    //
    //         const a_ = __a.to_lch();
    //         const b_ = __b.to_lch();
    //         const c1 = a_.c*a_.l*@as(f32, @floatFromInt(c1_));
    //         const c2 = b_.c*b_.l*@as(f32, @floatFromInt(c2_));
    //
    //         return
    //             a.partition_size > b.partition_size or
    //             c1 > c2;
    //     }
    // }.cmp);
    
    var primary = try gpa.alloc(Color.RGBA, means.len);
    defer gpa.free(primary);
    var complementary = try gpa.alloc(Color.RGBA, means.len);
    defer gpa.free(complementary);

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
    
    try pallet.print(gpa);


    // TODO: get input directory from args

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

    var cwd = std.fs.cwd();
    var in = try cwd.openDir("src/.ignore/template", .{});
    defer in.close();
    var out = try cwd.openDir("src/.ignore/out", .{});
    defer out.close();

    // try Parser.iterate_dir_generating_template(gpa, out, in, 0, pallet);







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
