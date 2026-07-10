# Deploying the reachability explorer

`index.html` in this folder is the interactive **reachability explorer** — a single, fully
self-contained HTML file. It embeds its own data (as an inline `const DATA = {…}` JSON payload) and
draws everything with inline JavaScript and SVG. **There is no build step and no server to run.**

## Fastest path: open it locally

You do not need to deploy anything to use it. Just open the file in any modern browser:

```bash
open deploy/index.html          # macOS
xdg-open deploy/index.html      # Linux
start deploy\index.html         # Windows
```

It works offline, from a `file://` URL — no internet connection required.

## Publish on GitHub Pages (free static hosting)

Because the explorer is one static file, GitHub Pages serves it as-is. Naming it `index.html` means the
folder's URL resolves straight to the explorer (no `/reachability_explorer.html` suffix needed).

1. **Push the repo** (this `deploy/` folder included) to GitHub, on the branch you want to publish from
   (usually `main`).
2. **Enable Pages.** In the repository on github.com, go to **Settings → Pages**.
   - Under **Build and deployment → Source**, choose **Deploy from a branch**.
   - Set **Branch** to `main` and the **folder** to **`/deploy`**, then **Save.**
     *(If your account only offers `/ (root)` and `/docs`, either move `index.html` to a top-level
     `docs/` folder and pick `/docs`, or publish from the repo root and pick `/ (root)`.)*
3. **Wait ~1–2 minutes** for the first build. The Pages panel then shows a green check and the live URL.

### Resulting URL

```
https://<your-username>.github.io/<your-repo-name>/
```

For example, if the GitHub user is `octocat` and the repo is `cell-state-reachability`:

```
https://octocat.github.io/cell-state-reachability/
```

Because the file is named `index.html`, that bare folder URL loads the explorer directly. Every push to
the publish branch re-deploys automatically.

## Why no build step

| Concern | Status |
|---|---|
| External JS/CSS/CDN dependencies | **None** — no `<script src>`, no `<link href>`, no `@import`. |
| Data fetch at runtime | **None** — the dataset is embedded inline as `const DATA`; no `fetch`/XHR. |
| Web fonts / external images | **None** — system fonts and inline SVG only. |
| Server / API needed | **No** — it is a static document. |

You can confirm self-containment yourself:

```bash
# should print nothing (no external references)
grep -Ei '<script[^>]*src=|<link[^>]*href=|https?://|@import|fetch\(' deploy/index.html
```

*(The only `url(...)` occurrences in the file are internal SVG fragment references like `url(#ah)` —
these point at `<defs>` inside the same document and are offline-safe.)*

## Custom domain (optional)

To serve it from your own domain, add a `CNAME` file next to `index.html` containing the domain (e.g.
`explorer.example.org`) and configure the DNS record per GitHub's
[custom-domain docs](https://docs.github.com/pages/configuring-a-custom-domain-for-your-github-pages-site).
