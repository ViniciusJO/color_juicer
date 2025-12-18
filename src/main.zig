const std = @import("std");

const stb_image = @import("stbi");
// const stb_image_write = @import("stbiw");
const stb_image_resize = @import("stbir");

const gpa = std.heap.page_allocator;

const Vec2 = [2]u8;
const Vec3 = [3]u8;
const Vec4 = [4]u8;

const Color = Vec4;
const HLV = struct { h: f32, s: f32, v: f32 };
const RGB = struct { r: f32, g: f32, b: f32 };
const Colors = std.ArrayList(Color);

const Mean = struct { color: Color, colors: Colors, dist: f64, partition_size: u64 = 0 };
const Means = std.ArrayList(Mean);

fn vec4_dist_sz(v1: Vec4, v2: [4]usize) f64 {
    const dx: f64 = @floatFromInt(v1[0] - v2[0]);
    const dy: f64 = @floatFromInt(v1[1] - v2[1]);
    const dz: f64 = @floatFromInt(v1[2] - v2[2]);
    const dw: f64 = @floatFromInt(v1[3] - v2[3]);
    return @sqrt(dx * dx + dy * dy + dz * dz + dw * dw);
}

fn vec4_dist(v1: Vec4, v2: Vec4) f64 {
    const dx: f64 = @floatFromInt(@abs(@as(i16, v1[0]) - @as(i16, v2[0])));
    const dy: f64 = @floatFromInt(@abs(@as(i16, v1[1]) - @as(i16, v2[1])));
    const dz: f64 = @floatFromInt(@abs(@as(i16, v1[2]) - @as(i16, v2[2])));
    const dw: f64 = @floatFromInt(@abs(@as(i16, v1[3]) - @as(i16, v2[3])));
    return @sqrt(dx * dx + dy * dy + dz * dz + dw * dw);
}

fn repartition(ms: *Means, cs: Colors) !*Means {
    for (ms.items) |*m| m.*.colors.clearRetainingCapacity(); //.clearRetainingCapacity();
    for (cs.items) |c| {
        var min = &ms.items[0];
        for (ms.items) |*m| {
            if (vec4_dist(c, m.*.color) < vec4_dist(c, min.color))
                min = m;
        }
        try min.colors.append(gpa, c);
        min.partition_size +|= 1;
    }
    return ms;
}

fn compute_means(ms: *Means, min_dist: f64) bool {
    var updated: usize = 0;
    for (ms.*.items) |*m| {
        if (m.*.dist <= min_dist) continue;

        var sum: [4]usize = .{0} ** 4;
        for (m.*.colors.items) |color| {
            sum[0] += color[0];
            sum[1] += color[1];
            sum[2] += color[2];
            sum[3] += color[3];
        }

        const color_res = if(m.*.colors.items.len > 0) Color {
            @intCast(sum[0] / m.*.colors.items.len),
            @intCast(sum[1] / m.*.colors.items.len),
            @intCast(sum[2] / m.*.colors.items.len),
            @intCast(sum[3] / m.*.colors.items.len)
        } else Color { 0 , 0, 0, 0 };

        m.*.dist = vec4_dist(m.*.color, color_res);
        m.*.color = color_res;
        updated += 1;
    }
    return updated == 0;
}

fn color_string(str: *const []const u8, c: Color) ![]u8 {
    return std.fmt.allocPrint(
        gpa,
        // "{d} + {d} = {d}",
        "\x1B[38;2;{d};{d};{d}m{s}\x1B[0m",
        .{ c[0], c[1], c[2], str.* },
    );
    // printf("\033[38;2;%d;%d;%dm%s\033[0m", r, g, b, str);
}

fn color_string_(str: *const []const u8, c: Color) ![]u8 {
    return std.fmt.allocPrint(
        gpa,
        // "{d} + {d} = {d}",
        "\x1B[38;2;{d};{d};{d}m{s}\x1B[0m",
        .{ c[0], c[1], c[2], str.* },
    );
    // printf("\033[38;2;%d;%d;%dm%s\033[0m", r, g, b, str);
}

pub fn writePPMImage(filename: []const u8, pixels: []const [3]u8) !void {
    if (pixels.len != 16) {
        return error.InvalidPixelCount;
    }

    const width = 4;
    const height = 4;
    const max_color = 255;

    var file = try std.fs.cwd().createFile(filename, .{});
    defer file.close();

    const writer = file.writer();

    // Write PPM header
    try writer.print("P3\n{} {}\n{}\n", .{ width, height, max_color });

    // Write pixel data (in row-major order)
    for (pixels) |pixel| {
        try writer.print("{} {} {}\n", .{ pixel[0], pixel[1], pixel[2] });
    }
}

fn downsampleImage(colors: Colors, width: u32, height: u32, factor: u32) !Colors {
    const newWidth = width / factor;
    const newHeight = height / factor;
    var downsampled = Colors.init(gpa);
    
    for (0..newHeight) |y| {
        for (0..newWidth) |x| {
            const origX = x * factor;
            const origY = y * factor;
            try downsampled.append(colors.items[origY * width + origX]);
        }
    }
    return downsampled;
}

fn rgbToHsv(r: f32, g: f32, b: f32) struct {h: f32, s: f32, v: f32} {
    const max = @max(r, @max(g, b));
    const min = @min(r, @min(g, b));
    const delta = max - min;

    var h: f32 = 0;
    if (delta != 0) {
        if (max == r) {
            h = 60 * @mod((g - b) / delta, 6);
        } else if (max == g) {
            h = 60 * ((b - r) / delta + 2);
        } else {
            h = 60 * ((r - g) / delta + 4);
        }
    }
    if (h < 0) h += 360;

    const s: f32 = if (max == 0) 0 else delta / max;
    const v: f32 = max;

    return .{ .h = h, .s = s, .v = v };
}

/// Convert HSV back to RGB.
/// Expects hue in [0,360), saturation and value in [0,1].
/// Returns tuple (r,g,b) each in [0,1].
fn hsvToRgb(h: f32, s: f32, v: f32) RGB {
    const c = v * s;
    const x = c * (1 - @abs(@mod(h / 60, 2) - 1));
    const m = v - c;

    var rp: RGB = .{ .r = 0, .g = 0, .b = 0 };

    if (h < 60) rp = .{ .r = c, .g = x, .b = 0 }
    else if (h < 120) rp = .{ .r = x, .g = c, .b = 0 }
    else if (h < 180) rp = .{ .r = 0, .g = c, .b = x }
    else if (h < 240) rp = .{ .r = 0, .g = x, .b = c }
    else if (h < 300) rp = .{ .r = x, .g = 0, .b = c }
    else rp = .{ .r = c, .g = 0, .b = x };

    return .{ .r = rp.r + m, .g = rp.g + m, .b = rp.b + m };
}

fn colorToRgb(c: Color) RGB {
    const rf = @as(f32, @floatFromInt(c[0])) / 255.0;
    const gf = @as(f32, @floatFromInt(c[1])) / 255.0;
    const bf = @as(f32, @floatFromInt(c[2])) / 255.0;
    return .{
        .r = rf,
        .g = gf,
        .b = bf,
    };
}

pub fn get_complementary(c: Color) Color {
    const rf = @as(f32, @floatFromInt(c[0])) / 255.0;
    const gf = @as(f32, @floatFromInt(c[1])) / 255.0;
    const bf = @as(f32, @floatFromInt(c[2])) / 255.0;

    const hsv = rgbToHsv(rf, gf, bf);

    // Add 180 degrees to hue and wrap around
    var h_comp = hsv.h + 180;
    if (h_comp >= 360) h_comp -= 360;

    const rgb_comp = hsvToRgb(h_comp, hsv.s, @abs(1.0-hsv.v));

    return .{
       @as(u8, @intFromFloat(std.math.clamp(rgb_comp.r * 255.0, 0, 255))),
       @as(u8, @intFromFloat(std.math.clamp(rgb_comp.g * 255.0, 0, 255))),
       @as(u8, @intFromFloat(std.math.clamp(rgb_comp.b * 255.0, 0, 255))),
       c[3], // Keep alpha unchanged
    };
}

fn uniformSeeds3D(ms: *Means) !*Means {
    const k = ms.items.len;

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

                ms.items[count].color[0] = r;
                ms.items[count].color[1] = g;
                ms.items[count].color[2] = b;
                count += 1;
            }
            if (count >= k) break;
        }
        if (count >= k) break;
    }

    return ms;
}

pub fn kmeanspp_init(m: *Means, pixels: *Colors) !*Means {
    const allocator = gpa;

    const n = pixels.items.len;
    const k = m.items.len;

    // Pick first centroid randomly
    const first_index = try randomIndex(n);
    const first_pixel = pixels.items[first_index];

    m.items[0].color = first_pixel;

    // Allocate distances array
    var distances = try allocator.alloc(f64, n);
    defer allocator.free(distances);

    // For each remaining centroid
    for (1..k) |i| {
        // For each pixel, find distance squared to nearest existing centroid
        for (pixels.items, 0..) |p, j| {
            var min_dist: f64 = std.math.inf(f64);
            for (m.items[0..i]) |mean| {
                const dist_sq = colorDistanceSq(mean.color, p);
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
        m.items[i].color = pixels.items[next_index];
    }

    return m;
}

fn colorDistanceSq(a: Color, b: Color) f64 {
    var sum: f64 = 0;
    for (a, b) |ai, bi| {
        const diff = @as(f64, @floatFromInt(ai)) - @as(f64, @floatFromInt(bi));
        sum += diff * diff;
    }
    return sum;
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

    errdefer std.debug.print("Usage {s}:\n\t{s} filepath <#means> <downsampling_factor>\nn", .{args[0], args[0]});

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

    // _ = res;
    // std.debug.print("res: {}, w: {d}, h: {d}, c: {d}\n", .{ img[0], w, h, c });
    var pixels = try std.ArrayList(Color).initCapacity(gpa, @as(usize, @intCast(w*h*c)));
        //.init(gpa);
    // _ = try pixels.addOne();
    var i: usize = 0;
    const size = w*h*c;
    while (i < size): (i += @intCast(c)) {
        try pixels.append(gpa, Color {
            img[i],
            img[i+1],
            img[i+2],
            255
        });
    }

    const new_w = @divTrunc(w, fac);
    const new_h = @divTrunc(h, fac);
    const new_c = 4;
    
    const downsampled = stb_image_resize.stbir_resize_uint8_linear(
        @ptrCast(pixels.items.ptr),
        w, h, w*4*@sizeOf(u8),
        null,
        new_w, new_h, new_w*4*@sizeOf(u8),
        new_c
    );

    if(null == downsampled) return error.ImageDownsamplingError;
    // else std.debug.print("Image Downsampled: ({}px, {}px, {} channels)\n", .{new_w, new_h, new_c});

    // _ = stb_image_write.stbi_write_png("downsampled.png", new_w, new_h, new_c, @ptrCast(downsampled), new_w*4*@sizeOf(u8));

    var pixels_1 = try Colors.initCapacity(gpa, 0);
    // _ = try pixels_1.addOne();
    var j: usize = 0;
    const size_1 = new_w*new_h*4;
    while (j < size_1): (j += @intCast(new_c)) {
        try pixels_1.append(gpa, Color {
            downsampled[j],
            downsampled[j+1],
            downsampled[j+2],
            255
        });
    }

    var means = try Means.initCapacity(gpa, 0);
    for (0..means_quant) |_| {
        try means.append(gpa, Mean {
            .color = Color { 0, 0, 0, 0xFF },
            .colors = try std.ArrayList(Color).initCapacity(gpa, 0),
            .dist = std.math.inf(f64),
            .partition_size = 0
        });
    }
    // _ = try uniformSeeds3D(&means);
    _ = uniformSeeds3D;
    _ = try kmeanspp_init(&means, &pixels_1);

    while (true) {
        _ = try repartition(&means, pixels_1);
        if (compute_means(&means, precision)) break;
    }

    for (means.items) |*m| {
        if(m.*.partition_size == 0) m.*.color = means.items[0].color;
    }
    

    // TODO: sort for proximity to background color (#000000FFF)
    std.sort.heap(Mean, means.items, {}, struct {
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

    var complementary = try Means.initCapacity(gpa, 0);
    for (means.items) |*m| {
        const comp = get_complementary(m.*.color);
        try complementary.append(gpa, Mean{
            .color = comp,
            .colors = undefined,
            .dist = 0,
            .partition_size = m.*.partition_size
        });
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
    for (complementary.items) |*m| {
        sum += m.*.partition_size;
        comp_mean[0] += m.*.color[0] * m.*.partition_size;
        comp_mean[1] += m.*.color[1] * m.*.partition_size;
        comp_mean[2] += m.*.color[2] * m.*.partition_size;
    }
    comp_mean[0] /= sum;
    comp_mean[1] /= sum;
    comp_mean[2] /= sum;

    const comp_m =  Color{ @intCast(comp_mean[0]), @intCast(comp_mean[1]), @intCast(comp_mean[2]), 255 };

    // for(0..4) |_| {
    //     std.debug.print("{s}", .{try color_string(&"████████████", comp_m)});
    //     std.debug.print("\n", .{});
    // }
    // std.debug.print("\n\n", .{});

    for (means.items) |*m| {
        var buff: [256]u8 = .{0}**256;
        const col = m.*.color;
        var writer = out_file.writer(&buff);
        var writeri = &writer.interface;
        try writeri.print("{} {} {}\n", .{col[0], col[1], col[2]});
        try writeri.flush();
    }
    
    const clrs = means.items;
    const cclrs = complementary.items;


    // try std.io.getStdOut().writer().print("prim: {s}\n", .{ try color_string_(&"██", clrs[0].color) });
    // try std.io.getStdOut().writer().print("sec: {s}\n", .{ try color_string_(&"██", clrs[1].color) });
    // try std.io.getStdOut().writer().print("cprim: {s}\n", .{ try color_string_(&"██", cclrs[0].color) });
    // try std.io.getStdOut().writer().print("csec: {s}\n", .{ try color_string_(&"██", cclrs[1].color) });
    // try std.io.getStdOut().writer().print("cont: {s}\n", .{ try color_string_(&"██", comp_m) });

    const color_dist_square = struct { fn color_dist(_clr: Color) u64 {
        const clr = [_]u64{ _clr[0], _clr[1], _clr[2], _clr[3] };
        return clr[0]*clr[0] + clr[1]*clr[1] + clr[2]*clr[2];
    }}.color_dist;

    std.sort.heap(Mean, means.items, {}, struct {
        pub fn cmp(_: void, _a: Mean, _b: Mean) bool {
            const a = colorToRgb(_a.color);
            const b = colorToRgb(_b.color);
            const c1_ = color_dist_square(_a.color);
            const c2_ = color_dist_square(_b.color);

            const a_ = rgbToHsv(a.r, a.g, a.b);
            const b_ = rgbToHsv(b.r, b.g, b.b);
            const c1 = a_.v*a_.s*@as(f32, @floatFromInt(c1_));
            const c2 = b_.v*b_.s*@as(f32, @floatFromInt(c2_));

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
        fn col(co: Color) [7]u8 {
            var ret: [7]u8 = undefined;
            var stream = std.io.fixedBufferStream(&ret);
            stream.writer().print("#{X:02}{X:02}{X:02}", .{ co[0], co[1], co[2] }) catch { ret = .{ '#', '0', '0', '0', '0', '0', '0' }; };
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
        .complementary  = color_to_rgb_str(comp_m)
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
        try out_writeri.print("pprim = #{X}{X}{X}\n", .{ clrs[0].color[0],  clrs[0].color[1], clrs[0].color[2] });
        try out_writeri.print("psec = #{X}{X}{X}\n", .{ clrs[1].color[0],  clrs[1].color[1], clrs[1].color[2] });
        try out_writeri.print("pterc = #{X}{X}{X}\n", .{ clrs[2].color[0],  clrs[2].color[1], clrs[2].color[2] });
        try out_writeri.print("pcont = #{X}{X}{X}\n", .{ comp_m[0], comp_m[1], comp_m[2] });

        try out_writeri.print("cprim = #{X}{X}{X}\n", .{ cclrs[0].color[0],  cclrs[0].color[1], cclrs[0].color[2] });
        try out_writeri.print("csec = #{X}{X}{X}\n", .{ cclrs[1].color[0],  cclrs[1].color[1], cclrs[1].color[2] });
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
    try stdout.print("pprim: {s} {s}\n", .{ try color_string_(&"██", clrs[0].color), color_to_rgb_str(clrs[0].color) });
    try stdout.print("psec:  {s} {s}\n", .{ try color_string_(&"██", clrs[1].color), color_to_rgb_str(clrs[1].color) });
    try stdout.print("pterc: {s} {s}\n", .{ try color_string_(&"██", clrs[2].color), color_to_rgb_str(clrs[2].color) });
    try stdout.print("pcont: {s} {s}\n", .{ try color_string_(&"██", comp_m), color_to_rgb_str(comp_m) });


    // for (means.items) |*m| {
    //     std.debug.print("{any}\n", .{m.*.color});
    //     std.debug.print("#{x}{x}{x}{x}\n", .{m.*.color[0], m.*.color[1], m.*.color[2], m.*.color[3]});
    //     // std.debug.print("\n{s}\n", .{try color_string(&"██", m.*.color)});
    //     try writer.print("{} {} {}\n", .{m.*.color[0], m.*.color[1], m.*.color[2]});
    // }
}
