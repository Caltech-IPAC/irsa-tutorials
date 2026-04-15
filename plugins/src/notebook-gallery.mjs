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


// MyST markdown templates — kept separate from logic for readability
const clickableCardMyst = (link, title, description) => [
  '```{card}',
  `:link: ${link}`,
  `:header: [${title} →](${link})`,
  description,
  '```',
].join('\n');

const unrecognizedCardMyst = (filePath) => [
  '```{card}',
  ':header: ⚠️ _Unrecognised notebook_',
  `Could not find \`${filePath}\``,
  '```',
].join('\n');

const metadataErrorMyst = (metadataPath) => [
  `:::{error} \`notebook-gallery\` error`,
  `Could not read metadata from \`${metadataPath}\``,
  ':::',
].join('\n');

const gridOfCardsMyst = (cards) => [
  '````{grid} 1 2 2 3',
  ...cards,
  '````',
].join('\n\n');
// End of MyST markdown templates


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
    // Resolve the metadata file path relative to the repo root and load it
    const metadataPath = path.resolve(REPO_ROOT, data.arg);
    const metadataByFile = loadMetadata(metadataPath);

    if (!metadataByFile) {
      return ctx.parseMyst(metadataErrorMyst(metadataPath)).children;
    }

    // Links in cards must be relative to the file that contains this directive
    const callerDir = path.dirname(path.resolve(vfile.path));

    // Parse the body: one notebook path per line, skip blank lines and comments
    const notebookFilePaths = data.body.split('\n')
      .map(l => l.trim())
      .filter(l => l && !l.startsWith('#'));

    // Build a MyST card string for each notebook path
    const cards = notebookFilePaths.map(filePath => {
      const meta = metadataByFile[filePath];
      if (!meta) return unrecognizedCardMyst(filePath);
      const link = path.relative(callerDir, path.resolve(REPO_ROOT, filePath));
      const title = meta.title ?? path.basename(filePath, path.extname(filePath)); // fall back to filename
      const description = meta.description ?? ''; // fall back to no description
      return clickableCardMyst(link, title, description);
    });

    // Wrap all cards in a grid, parse the combined MyST string, and return the
    // resulting AST nodes to be inserted in place of this directive.
    return ctx.parseMyst(gridOfCardsMyst(cards)).children;
  },
};

export default { name: 'IRSA Notebook Gallery', directives: [notebookGalleryDirective] };
