// Graphe des chapitres affiché en pyramide sur l'accueil : Terminale à la base, agrégation au sommet, un lien par
// prérequis. Les liens viennent de content/graphe.json (choix de structure, voir CLAUDE.md, « Format du contenu ») et des
// prérequis des cartes rédigées. Pur et testé : le composant Pyramide ne fait que dessiner la disposition calculée ici.

import { NIVEAUX, rangNiveau, type Chapitre, type Contenu, type Niveau } from './types'

/** Contenu de content/graphe.json : noms courts des chapitres et liens [prérequis, chapitre]. */
export interface FichierGraphe {
  courts: Record<string, string>
  liens: [string, string][]
}

/** Lien orienté : `de` est un prérequis de `vers`. */
export interface Lien {
  de: string
  vers: string
}

export const chapitresDe = (c: Contenu): Chapitre[] => c.themes.flatMap((t) => t.chapitres)

/** Liens de graphe.json, puis ceux qu'imposent les prérequis des cartes d'un chapitre à l'autre, sans doublon. */
export function liensChapitres(c: Contenu, g: FichierGraphe): Lien[] {
  const vus = new Set<string>()
  const sortie: Lien[] = []
  const ajouter = (de: string, vers: string) => {
    if (de === vers || vus.has(`${de}>${vers}`)) return
    vus.add(`${de}>${vers}`)
    sortie.push({ de, vers })
  }
  for (const [de, vers] of g.liens) ajouter(de, vers)
  const chapitreDe = new Map(c.cartes.map((k) => [k.id, k.chapitre]))
  for (const k of c.cartes)
    for (const p of k.prerequis) {
      const de = chapitreDe.get(p)
      if (de) ajouter(de, k.chapitre)
    }
  return sortie
}

/** Erreurs du fichier graphe.json et du graphe complet (vide si tout va bien). */
export function validerGraphe(c: Contenu, g: FichierGraphe): string[] {
  const err: string[] = []
  const niveaux = new Map(chapitresDe(c).map((ch) => [ch.id, ch.niveau]))
  for (const id of niveaux.keys()) if (!g.courts[id]?.trim()) err.push(`graphe.json : nom court manquant pour ${id}`)
  for (const id of Object.keys(g.courts)) if (!niveaux.has(id)) err.push(`graphe.json : chapitre inconnu ${id} (courts)`)
  const vus = new Set<string>()
  for (const lien of g.liens) {
    if (!Array.isArray(lien) || lien.length !== 2) {
      err.push(`graphe.json : lien mal formé ${JSON.stringify(lien)}`)
      continue
    }
    const [de, vers] = lien
    const nd = niveaux.get(de)
    const nv = niveaux.get(vers)
    if (!nd || !nv) err.push(`graphe.json : lien ${de} → ${vers} vers un chapitre inconnu`)
    else if (rangNiveau(nd) > rangNiveau(nv)) err.push(`graphe.json : lien ${de} → ${vers} qui descend de niveau`)
    if (de === vers) err.push(`graphe.json : lien de ${de} vers lui-même`)
    if (vus.has(`${de}>${vers}`)) err.push(`graphe.json : lien ${de} → ${vers} en double`)
    vus.add(`${de}>${vers}`)
  }
  const liens = liensChapitres(c, g).filter((l) => niveaux.has(l.de) && niveaux.has(l.vers))
  // Chaque chapitre est relié, et aucun cycle (un chapitre ne peut pas être son propre prérequis, même de loin).
  const relies = new Set(liens.flatMap((l) => [l.de, l.vers]))
  for (const id of niveaux.keys()) if (!relies.has(id)) err.push(`graphe : chapitre ${id} relié à rien`)
  const suivants = new Map<string, string[]>()
  for (const l of liens) suivants.set(l.de, [...(suivants.get(l.de) ?? []), l.vers])
  const etat = new Map<string, 'en-cours' | 'fini'>()
  const visiter = (id: string, chemin: string[]): void => {
    if (etat.get(id) === 'fini') return
    if (etat.get(id) === 'en-cours') {
      err.push(`graphe : cycle ${[...chemin.slice(chemin.indexOf(id)), id].join(' → ')}`)
      return
    }
    etat.set(id, 'en-cours')
    for (const s of suivants.get(id) ?? []) visiter(s, [...chemin, id])
    etat.set(id, 'fini')
  }
  for (const id of niveaux.keys()) visiter(id, [])
  return err
}

// Disposition, en unités du viewBox SVG.
export const ECART = 70
const INTERLIGNE = 64
const ENTRE_NIVEAUX = 28
const MARGE_HAUT = 36
const MARGE_GAUCHE = 76
const MARGE_DROITE = 24
/** Places d'une rangée : la pyramide s'élargit de haut en bas. */
const places = (rangee: number) => 4 + rangee

export interface NoeudPlace {
  id: string
  niveau: Niveau
  x: number
  y: number
}

/** Tranche d'un niveau : trapèze de yHaut à yBas, de demi-largeur demiHaut puis demiBas, centré sur la pyramide. */
export interface Bande {
  niveau: Niveau
  yHaut: number
  yBas: number
  demiHaut: number
  demiBas: number
}

export interface Pyramide {
  noeuds: NoeudPlace[]
  bandes: Bande[]
  largeur: number
  hauteur: number
  centre: number
}

/** Répartit n nœuds sur des rangées de capacités données, au prorata (plus forts restes). */
function repartir(n: number, capacites: number[]): number[] {
  const total = capacites.reduce((a, b) => a + b, 0)
  const parts = capacites.map((c) => (n * c) / total)
  const nombres = parts.map(Math.floor)
  const ordre = parts.map((p, i) => ({ i, reste: p - Math.floor(p) })).sort((a, b) => b.reste - a.reste || a.i - b.i)
  const manquants = n - nombres.reduce((a, b) => a + b, 0)
  for (let k = 0; k < manquants; k++) nombres[ordre[k].i]++
  return nombres
}

/**
 * Place chaque chapitre : un niveau par tranche (agrégation en haut), plusieurs rangées par tranche quand il le faut, puis
 * quelques passes de barycentres pour rapprocher chaque chapitre de ceux auxquels il est relié (moins de croisements).
 */
export function disposerPyramide(c: Contenu, liens: Lien[]): Pyramide {
  const chapitres = chapitresDe(c)
  const tranches = [...NIVEAUX].reverse().map((niveau) => ({ niveau, ids: chapitres.filter((ch) => ch.niveau === niveau).map((ch) => ch.id) }))
  const derniere = tranches.reduce((r, t) => {
    let k = 1
    while (Array.from({ length: k }, (_, i) => places(r + i)).reduce((a, b) => a + b, 0) < t.ids.length) k++
    return r + k
  }, 0)
  const centre = MARGE_GAUCHE + (places(derniere - 1) * ECART) / 2

  // Emplacements de chaque tranche (rangées décalées d'un quart d'écart pour ne pas empiler les nœuds), triés de gauche à droite.
  type Emplacement = { x: number; y: number }
  const emplacements: Emplacement[][] = []
  const bandes: Bande[] = []
  let rangee = 0
  tranches.forEach((t, it) => {
    let k = 1
    while (Array.from({ length: k }, (_, i) => places(rangee + i)).reduce((a, b) => a + b, 0) < t.ids.length) k++
    const capacites = Array.from({ length: k }, (_, i) => places(rangee + i))
    const nombres = repartir(t.ids.length, capacites)
    const y = (r: number) => MARGE_HAUT + r * INTERLIGNE + it * ENTRE_NIVEAUX
    const liste: Emplacement[] = []
    nombres.forEach((m, i) => {
      const decalage = k > 1 ? (i % 2 ? 1 : -1) * (ECART / 4) : 0
      for (let j = 0; j < m; j++) liste.push({ x: centre + (j - (m - 1) / 2) * ECART + decalage, y: y(rangee + i) })
    })
    emplacements.push(liste.sort((a, b) => a.x - b.x || a.y - b.y))
    bandes.push({
      niveau: t.niveau,
      yHaut: y(rangee) - INTERLIGNE / 2 + 4,
      yBas: y(rangee + k - 1) + INTERLIGNE / 2 + 4,
      demiHaut: (places(rangee - 0.5) * ECART) / 2,
      demiBas: (places(rangee + k - 0.5) * ECART) / 2,
    })
    rangee += k
  })

  // Ordre initial : celui de themes.json ; puis barycentres des voisins, en remontant puis en redescendant.
  const ordres = tranches.map((t) => [...t.ids])
  const voisins = new Map<string, string[]>()
  for (const l of liens) {
    voisins.set(l.de, [...(voisins.get(l.de) ?? []), l.vers])
    voisins.set(l.vers, [...(voisins.get(l.vers) ?? []), l.de])
  }
  const position = new Map<string, number>()
  const placer = (it: number) => ordres[it].forEach((id, j) => position.set(id, emplacements[it][j].x))
  ordres.forEach((_, it) => placer(it))
  const passe = (it: number) => {
    const cle = new Map(
      ordres[it].map((id) => {
        const xs = (voisins.get(id) ?? []).map((v) => position.get(v)).filter((x) => x !== undefined)
        return [id, xs.length ? xs.reduce((a, b) => a + b, 0) / xs.length : position.get(id)!]
      }),
    )
    ordres[it].sort((a, b) => cle.get(a)! - cle.get(b)!)
    placer(it)
  }
  const bas = ordres.map((_, i) => ordres.length - 1 - i)
  for (let n = 0; n < 4; n++) {
    bas.forEach(passe)
    ordres.forEach((_, it) => passe(it))
  }

  const noeuds = tranches.flatMap((t, it) => ordres[it].map((id, j) => ({ id, niveau: t.niveau, ...emplacements[it][j] })))
  const derniereBande = bandes[bandes.length - 1]
  return {
    noeuds,
    bandes,
    largeur: centre + derniereBande.demiBas + MARGE_DROITE,
    hauteur: derniereBande.yBas + 12,
    centre,
  }
}
