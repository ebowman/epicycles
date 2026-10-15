import { defineConfig, type Plugin } from "vite";
import { viteSingleFile } from "vite-plugin-singlefile";

// Drop the woff/ttf fallbacks from KaTeX's CSS so only woff2 fonts get inlined.
const katexWoff2Only: Plugin = {
  name: "katex-woff2-only",
  enforce: "pre",
  transform(code, id) {
    if (!id.includes("katex") || !id.endsWith(".css")) return null;
    return code.replace(/,url\([^)]*\.(?:woff|ttf)\) format\("(?:woff|truetype)"\)/g, "");
  },
};

// Single self-contained index.html (JS, CSS, KaTeX woff2 fonts, data all inlined).
export default defineConfig({
  base: "./",
  publicDir: false,
  plugins: [katexWoff2Only, viteSingleFile()],
  build: { outDir: "dist-single", emptyOutDir: true },
});
