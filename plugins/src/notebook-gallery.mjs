import fs from 'node:fs';
import path from 'node:path';
import process from 'node:process';
import YAML from 'yaml';

// myst build always runs from the repo root
const REPO_ROOT = process.cwd();

// Create a path-keyed cache so multiple directives calls sharing the same
// metadata file only trigger one disk read per build.
const metadataCache = new Map();

function loadMetadata(metadataPath) {
  if (metadataCache.has(metadataPath)) return metadataCache.get(metadataPath);
  try {
    const metadataText = fs.readFileSync(metadataPath, 'utf8');
    const metadataList = /\.json$/i.test(metadataPath)
      ? JSON.parse(metadataText) : YAML.parse(metadataText);
    // convert [{file, title, description}] to {file: {file, title, description}} for easy lookup
    const metadataMap = Object.fromEntries(metadataList.map(e => [e.file, e])); 
    metadataCache.set(metadataPath, metadataMap);
    return metadataMap;
  } catch {
    metadataCache.set(metadataPath, null);
    return null;
  }
}


/**
 * MyST directive: notebook-gallery
 *
 * Renders a list of notebooks as gallery: responsive grid of clickable cards.
 *
 * Arg:
 *   Path to a YAML or JSON metadata file, relative to the repo root.
 *   Must be a list of objects; each object must have at least the following fields:
 *     - title       {string}  Display title for the card header
 *     - file        {string}  Path to the notebook relative to the repo root
 *     - description {string}  One-sentence summary shown in the card body
 *
 * Body:
 *   One notebook file path (relative to the repo root) per line. This must 
 *   match the `file` field in the metadata. Lines starting with # are 
 *   treated as comments.
 *
 * Example:
 *   ```{notebook-gallery} metadata.yml
 *   tutorials/abc.md
 *   tutorials/def.md
 *   ```
 */
const notebookGalleryDirective = {
  name: 'notebook-gallery',
  arg: { type: String, required: true },
  body: { type: String, required: true },

  run(data, vfile, ctx) {
    const metadataPath = path.resolve(REPO_ROOT, data.arg);
    const metadataByFile = loadMetadata(metadataPath);

    if (!metadataByFile) {
      return ctx.parseMyst(`:::danger\nnotebook-gallery: could not read \`${metadataPath}\`\n:::`).children;
    }

    const callerDir = path.dirname(path.resolve(vfile.path));

    const cards = data.body.split('\n')
      .map(l => l.trim())
      .filter(l => l && !l.startsWith('#'))
      .map(filePath => {
        const meta = metadataByFile[filePath];
        if (!meta) return `:::warning\nnotebook-gallery: not found in metadata: \`${filePath}\`\n:::`;
        const link = path.relative(callerDir, path.resolve(REPO_ROOT, filePath));
        return `\`\`\`{card}\n:link: ${link}\n:header: [${meta.title} →](${link})\n${meta.description}\n\`\`\``;
      });

    const myst = `\`\`\`\`{grid} 1 2 2 3\n${cards.join('\n\n')}\n\`\`\`\``;
    return ctx.parseMyst(myst).children;
  },
};

export default { name: 'IRSA Notebook Gallery', directives: [notebookGalleryDirective] };
