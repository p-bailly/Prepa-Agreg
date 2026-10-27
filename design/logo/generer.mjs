// Logo de l'appli : « Le crayon et le carré de fin de démonstration », dessiné au crayon comme les logos de la famille
// Nimbus (Scrabble, Mes applis…) : papier #F6F1E6, graphite #2E3644, chaque trait repassé deux fois avec un léger tremblé,
// un aplat léger puis des hachures appuyées. Prépa Agreg prend le crayon rouge (le rouge des cours, #CE2029).
// Sans filtre : uniquement des tracés. Repris de design/logo/generer.mjs de p-bailly/Scrabble (mêmes outils de dessin).
//
// Écrit, sur la grille de l'icône adaptative Android (108 × 108, dessin dans le cercle sûr de 66) :
//  - design/logo/logo.svg (fond compris), repris par png.cjs pour les PNG ;
//  - les calques de l'icône adaptative Android (ic_launcher_foreground.xml, ic_launcher_monochrome.xml, fond papier)
//    et la petite icône des notifications (ic_stat_notify.xml, silhouette du carré et du crayon).
// Le tremblé vient d'un tirage pseudo-aléatoire à graine fixe : le dessin est toujours le même.
// Lancer : node design/logo/generer.mjs (puis node design/logo/png.cjs). Ne jamais éditer les fichiers produits à la main.

import { writeFileSync } from 'node:fs'

const RACINE = new URL('../../', import.meta.url)
const RES = new URL('android/app/src/main/res/', RACINE)

export const PAPIER = '#F6F1E6'
const GRAPHITE = '#2E3644'
const BLANC = '#FFFDF8'
const ROUGE_CLAIR = '#EE9A9F'
const ROUGE_APPUYE = '#CE2029'
const BOIS = '#E9CFA4'

// ---------- géométrie ----------
const rad = (d) => (d * Math.PI) / 180
const tourner = (pts, deg, cx, cy) => {
  const c = Math.cos(rad(deg)), s = Math.sin(rad(deg))
  return pts.map(([x, y]) => [cx + (x - cx) * c - (y - cy) * s, cy + (x - cx) * s + (y - cy) * c])
}
const decaler = (pts, dx, dy) => pts.map(([x, y]) => [x + dx, y + dy])

/** Ajoute des points le long de chaque côté (le tremblé est ainsi réparti sur tout le contour). */
function densifier(pts, ferme, pas = 1.3) {
  const sortie = []
  const n = ferme ? pts.length : pts.length - 1
  for (let i = 0; i < n; i++) {
    const a = pts[i], b = pts[(i + 1) % pts.length]
    const k = Math.max(1, Math.round(Math.hypot(b[0] - a[0], b[1] - a[1]) / pas))
    for (let j = 0; j < k; j++) sortie.push([a[0] + ((b[0] - a[0]) * j) / k, a[1] + ((b[1] - a[1]) * j) / k])
  }
  if (!ferme) sortie.push(pts[pts.length - 1])
  return sortie
}

function carre(cx, cy, cote, degres) {
  const h = cote / 2
  return tourner([[cx - h, cy - h], [cx + h, cy - h], [cx + h, cy + h], [cx - h, cy + h]], degres, cx, cy)
}

/** Point à l'abscisse `s` le long de l'axe du crayon, décalé de `t` perpendiculairement. */
function surAxe(debut, fin, s, t) {
  const L = Math.hypot(fin[0] - debut[0], fin[1] - debut[1])
  const u = [(fin[0] - debut[0]) / L, (fin[1] - debut[1]) / L], n = [-u[1], u[0]]
  return [debut[0] + u[0] * s + n[0] * t, debut[1] + u[1] * s + n[1] * t]
}

// ---------- le dessin, consigné en opérations puis écrit en SVG ou en vecteur Android ----------
const f1 = (n) => n.toFixed(1)
const D = (pts, ferme) => 'M' + pts.map((p) => `${f1(p[0])} ${f1(p[1])}`).join(' L') + (ferme ? ' Z' : '')

function melange(graine) {
  let a = graine >>> 0
  return () => {
    a |= 0
    a = (a + 0x6d2b79f5) | 0
    let t = Math.imul(a ^ (a >>> 15), 1 | a)
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296
  }
}

class Dessin {
  constructor(graine) {
    this.alea = melange(graine)
    this.ops = []
  }
  tremble(pts, a = 0.3) {
    return pts.map(([x, y]) => [x + (this.alea() * 2 - 1) * a, y + (this.alea() * 2 - 1) * a])
  }
  aplat(pts, couleur, alpha = 1) {
    this.ops.push({ k: 'aplat', d: D(pts, true), couleur, alpha })
  }
  /** Contour repassé deux fois : 100 % puis 55 % d'opacité, à 70 % de l'épaisseur. */
  trait(pts, { ferme = false, w = 2.6, couleur = GRAPHITE, j = 0.3 } = {}) {
    const dense = densifier(pts, ferme)
    for (const [largeur, alpha] of [[w, 1], [w * 0.7, 0.55]]) {
      this.ops.push({ k: 'trait', d: D(this.tremble(dense, j), ferme), couleur, w: largeur, alpha })
    }
  }
  contour(pts, o = {}) {
    this.trait(pts, { ferme: true, ...o })
  }
  /** Hachures appuyées, à −38° tous les 3,6, tenues dans la forme. */
  hachure(forme, couleur, { angle = -38, ecart = 3.6, w = 1.2, alpha = 0.85 } = {}) {
    const xs = forme.map((p) => p[0]), ys = forme.map((p) => p[1])
    const cx = (Math.min(...xs) + Math.max(...xs)) / 2, cy = (Math.min(...ys) + Math.max(...ys)) / 2
    const demi = Math.hypot(Math.max(...xs) - Math.min(...xs), Math.max(...ys) - Math.min(...ys)) / 2 + 3
    const a = rad(angle), u = [Math.cos(a), Math.sin(a)], n = [-Math.sin(a), Math.cos(a)]
    let d = ''
    for (let o = -demi; o <= demi; o += ecart + (this.alea() * 2 - 1) * 0.35) {
      const pente = (this.alea() * 2 - 1) * 0.12
      const l = demi + (this.alea() * 2 - 1) * 1.5
      const p0 = [cx + n[0] * o - u[0] * l + (this.alea() * 2 - 1), cy + n[1] * o - u[1] * l + (this.alea() * 2 - 1)]
      const p1 = [cx + n[0] * (o + pente * 10) + u[0] * l + (this.alea() * 2 - 1), cy + n[1] * (o + pente * 10) + u[1] * l + (this.alea() * 2 - 1)]
      d += `M${f1(p0[0])} ${f1(p0[1])} L${f1(p1[0])} ${f1(p1[1])} `
    }
    this.ops.push({ k: 'hachure', clip: D(forme, true), d: d.trim(), couleur, w, alpha })
  }
  ombre(forme) {
    this.hachure(forme, GRAPHITE, { ecart: 2.6, w: 0.9, alpha: 0.35 })
  }
}

// Le carré ∎ (fin de démonstration), plein, hachuré de rouge, légèrement penché.
const CARRE = { cx: 43.5, cy: 64.5, cote: 27, degres: -6 }
// Le crayon, de l'extrémité arrière (en haut à droite) à la mine, qui vient de finir le carré.
const ARRIERE = [80.5, 25.5], MINE = [58.5, 48.5]
const LARGEUR = 11

const dessin = new Dessin(20261027)
{
  const face = carre(CARRE.cx, CARRE.cy, CARRE.cote, CARRE.degres)
  dessin.ombre(decaler(face, 5, 6))
  dessin.aplat(face, ROUGE_CLAIR, 0.7)
  dessin.hachure(face, ROUGE_APPUYE, { ecart: 2.4, w: 1.4, alpha: 0.95 })
  dessin.contour(face, { w: 2.6 })

  const L = Math.hypot(MINE[0] - ARRIERE[0], MINE[1] - ARRIERE[1])
  const h = LARGEUR / 2
  const corpsFin = L - 8, boisFin = L - 2.4
  const p = (s, t) => surAxe(ARRIERE, MINE, s, t)
  const corps = [p(0, -h), p(corpsFin, -h), p(corpsFin, h), p(0, h)]
  const bois = [p(corpsFin, -h), p(boisFin, -h * 0.29), p(boisFin, h * 0.29), p(corpsFin, h)]
  const mine = [p(boisFin, -h * 0.29), p(L, 0), p(boisFin, h * 0.29)]
  const bague = [p(2.5, -h), p(5.5, -h), p(5.5, h), p(2.5, h)]
  dessin.ombre(decaler(corps, 4, 5))
  dessin.aplat(corps, ROUGE_CLAIR, 0.55)
  // Hachures dans le sens du crayon (de l'arrière vers la mine : 135°).
  dessin.hachure(corps, ROUGE_APPUYE, { angle: 135, ecart: 2.6, w: 1.2 })
  dessin.contour(corps, { w: 2.4 })
  dessin.aplat(bague, BLANC)
  dessin.contour(bague, { w: 1.6 })
  dessin.aplat(bois, BOIS)
  dessin.contour(bois, { w: 2.2 })
  dessin.aplat(mine, GRAPHITE)
}

// SVG : le même assemblage que les logos Nimbus (papier, aplats, hachures coupées par un clipPath, contours).
function corpsSvg(ops) {
  const defs = []
  const corps = ops.map((op) => {
    if (op.k === 'aplat') return `<path d="${op.d}" fill="${op.couleur}"${op.alpha < 1 ? ` fill-opacity="${op.alpha}"` : ''}></path>`
    if (op.k === 'trait') return `<path d="${op.d}" fill="none" stroke="${op.couleur}" stroke-opacity="${op.alpha.toFixed(2)}" stroke-width="${op.w.toFixed(2)}" stroke-linecap="round" stroke-linejoin="round"></path>`
    const id = `h${defs.length + 1}`
    defs.push(`<clipPath id="${id}"><path d="${op.clip}"></path></clipPath>`)
    return `<path d="${op.d}" stroke="${op.couleur}" stroke-width="${op.w}" stroke-linecap="round" stroke-opacity="${op.alpha}" clip-path="url(#${id})"></path>`
  })
  return `<defs>${defs.join('')}</defs>` + corps.join('')
}
const ETIQUETTE = 'Prépa Agreg : un crayon rouge qui vient de tracer le carré de fin de démonstration'
const svg = (viewBox, fond) =>
  `<svg xmlns="http://www.w3.org/2000/svg" viewBox="${viewBox}" role="img" aria-label="${ETIQUETTE}">${fond}${corpsSvg(dessin.ops)}</svg>\n`

writeFileSync(new URL('design/logo/logo.svg', RACINE), svg('0 0 108 108', `<rect width="108" height="108" fill="${PAPIER}"></rect>`))

// Vecteur Android : mêmes tracés, hachures coupées par un clip-path ; le calque monochrome garde traits et hachures.
const vecteur = (corps, taille = 108) => `<?xml version="1.0" encoding="utf-8"?>
<!-- Généré par design/logo/generer.mjs : ne pas modifier à la main. -->
<vector xmlns:android="http://schemas.android.com/apk/res/android"
    android:width="${taille}dp"
    android:height="${taille}dp"
    android:viewportWidth="108"
    android:viewportHeight="108">
${corps}
</vector>
`
function corpsAndroid(ops, mono) {
  const trait = (op) => `android:strokeColor="${mono ? '#FFFFFFFF' : op.couleur}" android:strokeAlpha="${op.alpha}" android:strokeWidth="${op.w}" android:strokeLineCap="round" android:strokeLineJoin="round"`
  return ops
    .map((op) => {
      if (op.k === 'aplat') return mono ? '' : `    <path android:fillColor="${op.couleur}"${op.alpha < 1 ? ` android:fillAlpha="${op.alpha}"` : ''} android:pathData="${op.d}" />\n`
      if (op.k === 'trait') return `    <path ${trait(op)} android:pathData="${op.d}" />\n`
      return `    <group>\n        <clip-path android:pathData="${op.clip}" />\n        <path ${trait(op)} android:pathData="${op.d}" />\n    </group>\n`
    })
    .join('')
    .trimEnd()
}
writeFileSync(new URL('drawable/ic_launcher_foreground.xml', RES), vecteur(corpsAndroid(dessin.ops, false)))
writeFileSync(new URL('drawable/ic_launcher_monochrome.xml', RES), vecteur(corpsAndroid(dessin.ops, true)))
writeFileSync(
  new URL('values/ic_launcher_background.xml', RES),
  `<?xml version="1.0" encoding="utf-8"?>\n<resources>\n    <color name="ic_launcher_background">${PAPIER}</color>\n</resources>\n`,
)

// Notification : silhouette d'une seule couleur, donc le carré plein et le crayon plein (sans tremblé), agrandis.
{
  const P = (pts) => D(pts, true)
  const L = Math.hypot(MINE[0] - ARRIERE[0], MINE[1] - ARRIERE[1])
  const h = LARGEUR / 2
  const p = (s, t) => surAxe(ARRIERE, MINE, s, t)
  const crayon = [p(0, -h), p(L - 8, -h), p(L, 0), p(L - 8, h), p(0, h)]
  const plein = `        <path android:fillColor="#FFFFFFFF" android:pathData="${P(carre(CARRE.cx, CARRE.cy, CARRE.cote, CARRE.degres))}" />\n` +
    `        <path android:fillColor="#FFFFFFFF" android:pathData="${P(crayon)}" />\n`
  writeFileSync(
    new URL('drawable/ic_stat_notify.xml', RES),
    vecteur(`    <group android:scaleX="1.45" android:scaleY="1.45" android:pivotX="54" android:pivotY="54">\n${plein}    </group>`, 24),
  )
}
console.log('logo.svg et icônes Android écrits')
