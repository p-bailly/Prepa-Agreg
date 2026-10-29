// Découpage du texte du contenu : paragraphes (ligne vide), listes (« - »), gras (**…**), formules $…$ et $$…$$,
// trous {{n}} des énoncés à trous. Pur et testé : le composant Texte ne fait que l'afficher, et le test du contenu
// s'en sert pour passer chaque formule à KaTeX. Voir CLAUDE.md, « Texte et formules ».

export type Morceau =
  | { k: 'texte'; t: string; gras?: boolean }
  | { k: 'math'; t: string; bloc: boolean }
  | { k: 'trou'; n: number }

export type Paragraphe = { k: 'p'; morceaux: Morceau[] } | { k: 'liste'; items: Morceau[][] } | { k: 'bloc'; t: string }

/** Découpe une ligne en texte, formules en ligne et trous. Une formule non refermée est gardée comme texte. */
export function decouperLigne(ligne: string): Morceau[] {
  const sortie: Morceau[] = []
  let i = 0
  let tampon = ''
  let gras = false
  const vider = () => {
    if (tampon) sortie.push(gras ? { k: 'texte', t: tampon, gras: true } : { k: 'texte', t: tampon })
    tampon = ''
  }
  while (i < ligne.length) {
    const c = ligne[i]
    if (c === '\\' && ligne[i + 1] === '$') {
      tampon += '$'
      i += 2
    } else if (c === '$') {
      const bloc = ligne[i + 1] === '$'
      const debut = i + (bloc ? 2 : 1)
      const fin = trouverFinMath(ligne, debut, bloc)
      if (fin < 0) {
        tampon += ligne.slice(i)
        break
      }
      vider()
      sortie.push({ k: 'math', t: ligne.slice(debut, fin).trim(), bloc })
      i = fin + (bloc ? 2 : 1)
    } else if (c === '*' && ligne[i + 1] === '*') {
      vider()
      gras = !gras
      i += 2
    } else if (c === '{' && ligne[i + 1] === '{') {
      const m = /^\{\{(\d+)\}\}/.exec(ligne.slice(i))
      if (m) {
        vider()
        sortie.push({ k: 'trou', n: Number(m[1]) })
        i += m[0].length
      } else {
        tampon += c
        i++
      }
    } else {
      tampon += c
      i++
    }
  }
  vider()
  return sortie
}

function trouverFinMath(ligne: string, debut: number, bloc: boolean): number {
  for (let j = debut; j < ligne.length; j++) {
    if (ligne[j] === '\\') {
      j++
      continue
    }
    if (ligne[j] === '$') {
      if (!bloc) return j
      if (ligne[j + 1] === '$') return j
    }
  }
  return -1
}

/**
 * Découpe un texte complet en paragraphes, listes et formules centrées. Dans un paragraphe, une suite de lignes « - »
 * forme une liste même collée à une ligne de texte (« Soit : » puis les cas) ; le texte autour reste en paragraphes.
 */
export function decouperTexte(texte: string): Paragraphe[] {
  const sortie: Paragraphe[] = []
  const estItem = (l: string) => /^\s*- /.test(l)
  for (const brut of texte.split(/\n\s*\n/)) {
    const lignes = brut.split('\n').map((l) => l.trimEnd()).filter((l) => l.trim() !== '')
    if (!lignes.length) continue
    // Une formule centrée seule sur ses lignes : bloc à part (elle défile dans son cadre si elle est trop large).
    const entier = lignes.join('\n').trim()
    if (/^\$\$[\s\S]*\$\$$/.test(entier) && entier.indexOf('$$', 2) === entier.length - 2) {
      sortie.push({ k: 'bloc', t: entier.slice(2, -2).trim() })
      continue
    }
    for (let i = 0; i < lignes.length; ) {
      let j = i
      while (j < lignes.length && estItem(lignes[j]) === estItem(lignes[i])) j++
      const groupe = lignes.slice(i, j)
      if (estItem(lignes[i])) sortie.push({ k: 'liste', items: groupe.map((l) => decouperLigne(l.replace(/^\s*- /, ''))) })
      else sortie.push({ k: 'p', morceaux: decouperLigne(groupe.join(' ')) })
      i = j
    }
  }
  return sortie
}

/** Toutes les formules d'un texte (pour le contrôle KaTeX). */
export function formulesDe(texte: string): { t: string; bloc: boolean }[] {
  const f: { t: string; bloc: boolean }[] = []
  for (const p of decouperTexte(texte)) {
    if (p.k === 'bloc') f.push({ t: p.t, bloc: true })
    const lignes = p.k === 'p' ? [p.morceaux] : p.k === 'liste' ? p.items : []
    for (const l of lignes) for (const m of l) if (m.k === 'math') f.push({ t: m.t, bloc: m.bloc })
  }
  return f
}

/** Numéros des trous {{n}} présents dans un texte, dans l'ordre d'apparition. */
export function trousDe(texte: string): number[] {
  return [...texte.matchAll(/\{\{(\d+)\}\}/g)].map((m) => Number(m[1]))
}

/** Numéros des trous vraiment affichés (hors formules) : un trou dans une formule ne serait pas reconnu. */
export function trousAffiches(texte: string): number[] {
  const n: number[] = []
  for (const p of decouperTexte(texte)) {
    const lignes = p.k === 'p' ? [p.morceaux] : p.k === 'liste' ? p.items : []
    for (const l of lignes) for (const m of l) if (m.k === 'trou') n.push(m.n)
  }
  return n
}
