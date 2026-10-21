# Presentation

Slides for the conference talk on the manuscript (PO-26-00422), built with
[reveal.js](https://revealjs.com) 6.0.2 ("full setup": the dev server gives
live reload while editing).

## Run locally

Requires Node.js >= 18.

```sh
cd presentation
npm install      # first time only
npm start        # serves http://localhost:8000 and reloads on save
```

## Where things go

```
index.html    the slides (one <section> per slide; nest sections for vertical slides)
assets/       figures and media used in slides, e.g. copies from ../output/figures
dist/         reveal.js runtime, themes and plugins (do not edit; regenerate with npm run build)
css/ js/ plugin/   reveal.js sources, only needed if you change the framework itself
```

Useful while editing:

- Theme: swap `dist/theme/white.css` in `index.html` for any file in `dist/theme/`.
- Markdown slides: `<section data-markdown>`, see the example slide in `index.html`.
- Math: KaTeX via the math plugin, write `\( ... \)` inline or `\[ ... \]` display.
- Speaker notes: `<aside class="notes">` inside a slide; press `S` during the talk.
- PDF export: open `http://localhost:8000/?print-pdf` in Chrome and print to PDF.
- Full docs: https://revealjs.com

## What is vendored here

A copy of the reveal.js repository with its own `.git`, `test/`, `react/`,
`examples/`, `demo.html` and CI config removed. `LICENSE` is reveal.js's MIT
license and applies to the framework files only. `node_modules/` is ignored by
git; `npm install` recreates it from `package-lock.json`.
