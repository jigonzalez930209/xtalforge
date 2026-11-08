import { defineConfig } from 'vitepress'

const roadmapItems = [
  { text: 'Overview', link: '/roadmap/' },
  { text: 'Phase 0 — Foundations', link: '/roadmap/phase-00-foundations' },
  { text: 'Phase 0 — Remaining work', link: '/roadmap/phase-00-status' },
  { text: 'Phase 1 — Core', link: '/roadmap/phase-01-core' },
  { text: 'Phase 1 — Methodology', link: '/roadmap/phase-01-methodology' },
  { text: 'Phase 2 — Symmetry', link: '/roadmap/phase-02-symmetry' },
  { text: 'Phase 3 — Formats, commands, view data', link: '/roadmap/phase-03-formats-commands' },
  { text: 'Phase 4 — Force fields', link: '/roadmap/phase-04-force-fields' },
  { text: 'Phase 5 — Analysis', link: '/roadmap/phase-05-analysis' },
  { text: 'Phase 6 — Builders', link: '/roadmap/phase-06-builders' },
  { text: 'Phase 7 — Powder refinement', link: '/roadmap/phase-07-powder' },
  { text: 'Phase 8 — Agent and distribution', link: '/roadmap/phase-08-distribution' },
]

export default defineConfig({
  title: 'xtalforge',
  description: 'Pure Rust library for building, analysing and exporting crystal structures',
  cleanUrls: true,
  lastUpdated: true,
  markdown: { math: true },
  editLink: {
    pattern: 'https://github.com/jigonzalez930209/xtalforge/edit/main/docs/:path',
    text: 'Edit this page on GitHub',
  },
  themeConfig: {
    nav: [
      { text: 'Guide', link: '/guide/' },
      { text: 'API', link: '/api/' },
      { text: 'Algorithms', link: '/algorithms/' },
      { text: 'Demos', link: '/demos/' },
      { text: 'Roadmap', link: '/roadmap/' },
      { text: 'Changelog', link: '/changelog' },
    ],
    sidebar: {
      '/guide/': [
        {
          text: 'Guide',
          items: [
            { text: 'Introduction', link: '/guide/' },
            { text: 'Contributing', link: '/guide/contributing' },
            { text: 'Cargo features', link: '/guide/features' },
            { text: 'Data bundles', link: '/guide/data' },
          ],
        },
      ],
      '/api/': [
        {
          text: 'API',
          items: [
            { text: 'Overview', link: '/api/' },
            { text: 'Binding architecture', link: '/api/bindings' },
            { text: 'Core API', link: '/api/core' },
            { text: 'Elements', link: '/api/elements' },
            { text: 'Hello-world API', link: '/api/hello-world' },
            { text: 'Surface parity matrix', link: '/api/parity' },
          ],
        },
      ],
      '/algorithms/': [
        {
          text: 'Algorithms',
          items: [
            { text: 'Overview', link: '/algorithms/' },
            { text: 'Validation and tolerances', link: '/algorithms/validation' },
            { text: 'Golden data', link: '/algorithms/golden-data' },
            { text: 'Symmetry backend: moyo', link: '/algorithms/symmetry-backend' },
            { text: 'sci-form reuse', link: '/algorithms/sci-form-reuse' },
            { text: 'Dependency replacement map', link: '/algorithms/dependency-map' },
          ],
        },
      ],
      '/demos/': [
        {
          text: 'Demos',
          items: [
            { text: 'Overview', link: '/demos/' },
            { text: 'Cell volume', link: '/demos/cell-volume' },
          ],
        },
      ],
      '/roadmap/': [{ text: 'Roadmap', items: roadmapItems }],
    },
    search: { provider: 'local' },
    socialLinks: [
      { icon: 'github', link: 'https://github.com/jigonzalez930209/xtalforge' },
    ],
  },
})
