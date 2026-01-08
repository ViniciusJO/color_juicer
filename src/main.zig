const std = @import("std");
const Color = @import("color.zig");
const Vec = @import("vec.zig");

const stb_image = @import("stbi");
// const stb_image_write = @import("stbiw");
const stb_image_resize = @import("stbir");

const gpa = std.heap.page_allocator;

const Mean = struct { color: Color.RGBA, colors: std.ArrayList(Color.RGBA), dist: f64, partition_size: u64 = 0 };

fn repartition(ms: *[]Mean, pixels: *[]const [4]u8) !*[]Mean {
    for (ms.*) |*m| m.*.colors.clearRetainingCapacity(); //.clearRetainingCapacity();
    for (pixels.*) |c| {
        var min = &(ms.*[0]);
        for (ms.*) |*m| {
            if (Vec.vec4_dist(c, m.color.to_vec4()) < Vec.vec4_dist(c, min.color.to_vec4())) min = m;
        }
        try min.colors.append(gpa, Color.RGBA.init(.{ .vec = c }));
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
            sum[0] += color.r;
            sum[1] += color.g;
            sum[2] += color.b;
            sum[3] += color.a;
        }

        const color_res = if(m.*.colors.items.len > 0) Vec.Vec4 {
            @intCast(sum[0] / m.*.colors.items.len),
            @intCast(sum[1] / m.*.colors.items.len),
            @intCast(sum[2] / m.*.colors.items.len),
            @intCast(sum[3] / m.*.colors.items.len)
        } else Vec.Vec4{ 0 , 0, 0, 0 };

        m.*.dist = Vec.vec4_dist(m.*.color.to_vec4(), color_res);
        m.*.color = Color.RGBA.init(.{ .vec = color_res });
        updated += 1;
    }
    return updated == 0;
}

fn uniformSeeds3D(ms: *[]Mean) !*[]Mean {
    const k = ms.len;

    const k_f: f64 = @floatFromInt(k);
    const divs_f = std.math.cbrt(k_f);
    const divs_unclamped: usize = @intFromFloat(@ceil(divs_f));
    const divs = std.math.clamp(divs_unclamped, 1, 256);

    const step = 256.0 / @as(f64, @floatFromInt(divs));

    var count: usize = 0;
    for (0..divs) |i| {
        for (0..divs) |j| {
            for (0..divs) |l| {
                if (count >= k) break;

                const r = @as(u8, @intFromFloat(std.math.clamp((@as(f64, @floatFromInt(i)) + 0.5) * step, 0, 255)));
                const g = @as(u8, @intFromFloat(std.math.clamp((@as(f64, @floatFromInt(j)) + 0.5) * step, 0, 255)));
                const b = @as(u8, @intFromFloat(std.math.clamp((@as(f64, @floatFromInt(l)) + 0.5) * step, 0, 255)));

                ms[count].color[0] = r;
                ms[count].color[1] = g;
                ms[count].color[2] = b;
                count += 1;
            }
            if (count >= k) break;
        }
        if (count >= k) break;
    }

    return ms;
}

pub fn kmeanspp_init(m: *[]Mean, pixels: *[]const [4]u8) !*[]Mean {
    const allocator = gpa;

    const n = pixels.len;
    const k = m.len;

    // Pick first centroid randomly
    const first_index = try randomIndex(n);
    const first_pixel = pixels.*[first_index];

    m.*[0].color = Color.RGBA.init(.{ .vec = first_pixel });

    // Allocate distances array
    var distances = try allocator.alloc(f64, n);
    defer allocator.free(distances);

    // For each remaining centroid
    for (1..k) |i| {
        // For each pixel, find distance squared to nearest existing centroid
        for (pixels.*, 0..) |p, j| {
            var min_dist: f64 = std.math.inf(f64);
            for (m.*[0..i]) |mean| {
                const dist_sq = colorDistanceSq(mean.color, Color.RGBA.init(.{ .vec = p }));
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
        m.*[i].color = Color.RGBA.init(.{ .vec = pixels.*[next_index] });
    }

    return m;
}

fn colorDistanceSq(a: Color.RGBA, b: Color.RGBA) f64 {
    const square = struct { pub fn sq(x: f64) f64 { return x*x; } }.sq;
    return
        square(@as(f64, @floatFromInt(a.r)) - @as(f64, @floatFromInt(b.r))) + 
        square(@as(f64, @floatFromInt(a.g)) - @as(f64, @floatFromInt(b.g))) + 
        square(@as(f64, @floatFromInt(a.b)) - @as(f64, @floatFromInt(b.b))) + 
        square(@as(f64, @floatFromInt(a.a)) - @as(f64, @floatFromInt(b.a)));
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

fn usage(args: [][:0]u8) void {
    std.debug.print("Usage {s}:\n\t{s} filepath <#means> <downsampling_factor>\nn", .{args[0], args[0]});
}

pub fn main() !void {
    const args = try std.process.argsAlloc(gpa);

    // const f = try std.fs.cwd().createFile("test.txt", .{});
    // defer f.close();
    //
    // var __b: [256]u8 = undefined;
    // var v = f.writer(&__b);
    // var in = &v.interface;
    //
    // try in.print("HEY\n", .{});
    // try in.flush();
    //
    // if(true) return;

    // try std.io.getStdOut().writer().print("\n{}: {s}\n", .{ args.len, args });

    errdefer usage(args);

    const img_path = if(args.len >= 1) args[1] else return error.NoFilepath;
    const means_quant = if(args.len >= 2) std.fmt.parseInt(usize, args[2], 10) catch return error.InvalidArgument else 10;
    const fac = if(args.len >= 3) std.fmt.parseInt(u8, args[3], 10) catch return error.InvalidArgument else 4;
    const precision = 0.01;

    var w: c_int = undefined;
    var h: c_int = undefined;
    var c: c_int = undefined;
    const img = stb_image.stbi_load(img_path, &w, &h, &c, 0);
    if(null == img) return error.ImageLoadError;
    // else std.debug.print("\nImage Loaded: {s} ({}px, {}px, {} channels)\n", .{img_path, w, h, c});

    const new_w = @divTrunc(w, fac);
    const new_h = @divTrunc(h, fac);
    const new_c = 3;
    
    const downsampled = stb_image_resize.stbir_resize_uint8_linear(
        img,
        w, h, w*c*@sizeOf(u8),
        null,
        new_w, new_h, new_w*4*@sizeOf(u8),
        new_c
    );
    if(null == downsampled) return error.ImageDownsamplingError;

    const down_size: usize = @intCast(new_w*new_h*new_c);
    var down_slice: []const u8 = downsampled[0..down_size];

    // else std.debug.print("Image Downsampled: ({}px, {}px, {} channels)\n", .{new_w, new_h, new_c});

    // _ = stb_image_write.stbi_write_png("downsampled.png", new_w, new_h, new_c, @ptrCast(downsampled), new_w*4*@sizeOf(u8));

    // var pixels_1 = try Color.Colors.initCapacity(gpa, 0);
    // // _ = try pixels_1.addOne();
    // var j: usize = 0;
    // const size_1 = new_w*new_h*4;
    // while (j < size_1): (j += @intCast(new_c)) {
    //     try pixels_1.append(gpa, Color.ColorInt {
    //         downsampled[j],
    //         downsampled[j+1],
    //         downsampled[j+2],
    //         255
    //     });
    // }

    // var means = try []Mean.initCapacity(gpa, 0);
    var means = try gpa.alloc(Mean, means_quant);
    // _ = &means_slice;
    // for (0..means_quant) |_| {
    //     try means.append(gpa, Mean {
    //         .color = Color.ColorInt { 0, 0, 0, 0xFF },
    //         .colors = try std.ArrayList(Color.ColorInt).initCapacity(gpa, 0),
    //         .dist = std.math.inf(f64),
    //         .partition_size = 0
    //     });
    // }

    for(means) |*m| {
        m.color = Color.RGBA.init(.{ .hex = 0x000000FF });
        m.colors = try std.ArrayList(Color.RGBA).initCapacity(gpa, 0);
        m.dist = std.math.inf(f64);
        m.partition_size = 0;
    }
    // _ = try uniformSeeds3D(&means);
    _ = uniformSeeds3D;
    _ = try kmeanspp_init(&means, &down_slice);

    while (true) {
        _ = try repartition(&means, &down_slice);
        if (compute_means(&means, precision)) break;
    }

    for (means) |*m| {
        if(m.*.partition_size == 0) m.*.color = means[0].color;
    }
    

    // TODO: sort for proximity to background color (#000000FFF)
    std.sort.heap(Mean, means, {}, struct {
        pub fn cmp(_: void, a: Mean, b: Mean) bool {
            return a.partition_size > b.partition_size;
        }
    }.cmp);

    const home_dir_path = try std.process.getEnvVarOwned(gpa, "HOME");
    var home = try std.fs.openDirAbsolute(home_dir_path, .{});
    defer home.close();
    var cache = try home.openDir(".cache/", .{});
    defer cache.close();
    var out_file = try cache.createFile("colours", .{});
    defer out_file.close();
    var out_ini_file = try cache.createFile("dyn_colors.ini", .{});
    defer out_ini_file.close();
    var out_i3_file = try cache.createFile("i3_colors", .{});
    defer out_i3_file.close();

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

    // var complementary = try []Mean.initCapacity(gpa, 0);
    var complementary = try gpa.alloc(Mean, means.len);
    for (means, 0..) |*m, i| {
        complementary[i] = Mean{
            .color = m.color.complementary(),
            .colors = undefined,
            .dist = 0,
            .partition_size = m.*.partition_size
        };
    }

    // for(0..4) |_| {
    //     for (complementary.items) |*m| {
    //         const col = m.*.color;
    //         std.debug.print("{s}", .{try color_string(&"████████████", col)});
    //     }
    //     std.debug.print("\n", .{});
    // }
    // for (complementary.items) |*m| {
    //     std.debug.print("{d:^12}", .{ m.*.partition_size });
    // }
    // std.debug.print("\n\n", .{});

    var comp_mean: [3]u64 = .{0, 0, 0};
    var sum: u64 = 0;
    for (complementary) |*m| {
        sum += m.*.partition_size;
        comp_mean[0] += m.*.color.r * m.*.partition_size;
        comp_mean[1] += m.*.color.g * m.*.partition_size;
        comp_mean[2] += m.*.color.b * m.*.partition_size;
    }
    comp_mean[0] /= sum;
    comp_mean[1] /= sum;
    comp_mean[2] /= sum;

    const comp_m =  Vec.Vec4{ @intCast(comp_mean[0]), @intCast(comp_mean[1]), @intCast(comp_mean[2]), 255 };

    // for(0..4) |_| {
    //     std.debug.print("{s}", .{try color_string(&"████████████", comp_m)});
    //     std.debug.print("\n", .{});
    // }
    // std.debug.print("\n\n", .{});

    for (means) |*m| {
        var buff: [256]u8 = .{0}**256;
        const col = m.*.color;
        var writer = out_file.writer(&buff);
        var writeri = &writer.interface;
        try writeri.print("{} {} {}\n", .{col.r, col.g, col.b});
        try writeri.flush();
    }
    
    const clrs = means;
    const cclrs = complementary;


    // try std.io.getStdOut().writer().print("prim: {s}\n", .{ try color_string_(&"██", clrs[0].color) });
    // try std.io.getStdOut().writer().print("sec: {s}\n", .{ try color_string_(&"██", clrs[1].color) });
    // try std.io.getStdOut().writer().print("cprim: {s}\n", .{ try color_string_(&"██", cclrs[0].color) });
    // try std.io.getStdOut().writer().print("csec: {s}\n", .{ try color_string_(&"██", cclrs[1].color) });
    // try std.io.getStdOut().writer().print("cont: {s}\n", .{ try color_string_(&"██", comp_m) });

    const color_dist_square = struct { fn color_dist(_clr: Color.RGBA) u64 {
        const clr = [_]u64{ _clr.r, _clr.g, _clr.b, _clr.a };
        return clr[0]*clr[0] + clr[1]*clr[1] + clr[2]*clr[2];
    }}.color_dist;

    std.sort.heap(Mean, means, {}, struct {
        pub fn cmp(_: void, _a: Mean, _b: Mean) bool {
            const a = _a.color;
            const b = _b.color;
            const c1_ = color_dist_square(_a.color);
            const c2_ = color_dist_square(_b.color);

            const a_ = a.to_lch();
            const b_ = b.to_lch();
            const c1 = a_.c*a_.l*@as(f32, @floatFromInt(c1_));
            const c2 = b_.c*b_.l*@as(f32, @floatFromInt(c2_));

            return c1 > c2;
        }
    }.cmp);

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

    const color_to_rgb_str = struct {
        fn col(co: Color.RGBA) [7]u8 {
            var ret: [7]u8 = undefined;
            var stream = std.io.fixedBufferStream(&ret);
            stream.writer().print("#{X:02}{X:02}{X:02}", .{ co.r, co.g, co.b }) catch { ret = .{ '#', '0', '0', '0', '0', '0', '0' }; };
            return ret;
        }
    }.col;

    const Colorscheme = struct {
        primary: [7]u8,
        secondary: [7]u8,
        terciary: [7]u8,
        complementary: [7]u8
    };

    const cs = Colorscheme {
        .primary  = color_to_rgb_str(clrs[0].color),
        .secondary  = color_to_rgb_str(clrs[0].color),
        .terciary  = color_to_rgb_str(clrs[0].color),
        .complementary  = color_to_rgb_str(Color.RGBA.init(.{ .vec = comp_m })),
    };

    // borda | fundo título | texto título | indicador | texto título (estado inverso)
    {
        var out_i3_writer = out_i3_file.writer(&.{});
        var out_i3_writeri = &out_i3_writer.interface;
        try out_i3_writeri.print("\nclient.focused {s} {s} {s} {s} {s}\n", .{
            cs.complementary, // title border
            cs.complementary, // title background
            "#000000",        // title text
            cs.terciary     , // indicator
            cs.complementary  // border
        });
    }

    // var out_ini_file = try std.fs.cwd().createFile("dyn_cpçprs.ini", .{});
    // defer out_ini_file.close();
    // // {
    //     // var bu = [1]u8{0}**255;
    //     // var ws = std.Io.fixedBufferStream(&bu);
    //     // var out_writer = ws.writer();
    //
        var out_writer = out_ini_file.writer(&.{});
        var out_writeri = &out_writer.interface;
        try out_writeri.print("[dyn_colors]\n", .{});
        try out_writeri.print("pprim = #{X}{X}{X}\n", .{ clrs[0].color.r,  clrs[0].color.g, clrs[0].color.b });
        try out_writeri.print("psec = #{X}{X}{X}\n", .{ clrs[1].color.r,  clrs[1].color.g, clrs[1].color.b });
        try out_writeri.print("pterc = #{X}{X}{X}\n", .{ clrs[2].color.r,  clrs[2].color.g, clrs[2].color.b });
        try out_writeri.print("pcont = #{X}{X}{X}\n", .{ comp_m[0], comp_m[1], comp_m[2] });

        try out_writeri.print("cprim = #{X}{X}{X}\n", .{ cclrs[0].color.r,  cclrs[0].color.g, cclrs[0].color.b });
        try out_writeri.print("csec = #{X}{X}{X}\n", .{ cclrs[1].color.r,  cclrs[1].color.g, cclrs[1].color.b });
        // try out_writeri.print("cont = #{X}{X}{X}\n", .{ comp_m[0], comp_m[1], comp_m[2] });

        // try out_out_1.print("pprim = #{X}{X}{X}\n", .{ clrs[0].color[0],  clrs[0].color[1], clrs[0].color[2] });
        // try out_out_1.print("psec = #{X}{X}{X}\n", .{ clrs[1].color[0],  clrs[1].color[1], clrs[1].color[2] });
        // try out_out_1.print("pterc = #{X}{X}{X}\n", .{ clrs[2].color[0],  clrs[2].color[1], clrs[2].color[2] });
        // try out_out_1.print("pcont = #{X}{X}{X}\n", .{ comp_m[0], comp_m[1], comp_m[2] });

        // std.debug.print("{s}", .{ bu });
        // std.debug.print("{s}", .{ bu });

        // std.debug.print("[dyn_colors]\n", .{});
        // std.debug.print("prim = #{X}{X}{X}\n", .{ clrs[0].color[0],  clrs[0].color[1], clrs[0].color[2] });
        // std.debug.print("sec = #{X}{X}{X}\n", .{ clrs[1].color[0],  clrs[1].color[1], clrs[1].color[2] });
        // std.debug.print("cprim = #{X}{X}{X}\n", .{ cclrs[0].color[0],  cclrs[0].color[1], cclrs[0].color[2] });
        // std.debug.print("csec = #{X}{X}{X}\n", .{ cclrs[1].color[0],  cclrs[1].color[1], cclrs[1].color[2] });
        // std.debug.print("cont = #{X}{X}{X}\n", .{ comp_m[0], comp_m[1], comp_m[2] });
        // try out_writer.flush();
    // }

    var stdout_writer = std.fs.File.stdout().writer(&.{});
    const stdout = &stdout_writer.interface;
    // try stdout.print("pprim: {s} {s}\n", .{ try color_string_(&"██", clrs[0].color), color_to_rgb_str(clrs[0].color) });
    // try stdout.print("psec:  {s} {s}\n", .{ try color_string_(&"██", clrs[1].color), color_to_rgb_str(clrs[1].color) });
    // try stdout.print("pterc: {s} {s}\n", .{ try color_string_(&"██", clrs[2].color), color_to_rgb_str(clrs[2].color) });
    // try stdout.print("pcont: {s} {s}\n", .{ try color_string_(&"██", comp_m), color_to_rgb_str(comp_m) });

    try stdout.print("pprim: {s} {s}\n", .{ try clrs[0].color.colorizer(gpa, &"██"), color_to_rgb_str(clrs[0].color) });
    try stdout.print("psec:  {s} {s}\n", .{ try clrs[1].color.colorizer(gpa, &"██"), color_to_rgb_str(clrs[1].color) });
    try stdout.print("pterc: {s} {s}\n", .{ try clrs[2].color.colorizer(gpa, &"██"), color_to_rgb_str(clrs[2].color) });
    try stdout.print("pcont: {s} {s}\n", .{ try Color.RGBA.init(.{ .vec = comp_m }).colorizer(gpa, &"██"), color_to_rgb_str(Color.RGBA.init(.{ .vec = comp_m })) });


    // for (means.items) |*m| {
    //     std.debug.print("{any}\n", .{m.*.color});
    //     std.debug.print("#{x}{x}{x}{x}\n", .{m.*.color[0], m.*.color[1], m.*.color[2], m.*.color[3]});
    //     // std.debug.print("\n{s}\n", .{try color_string(&"██", m.*.color)});
    //     try writer.print("{} {} {}\n", .{m.*.color[0], m.*.color[1], m.*.color[2]});
    // }
}
