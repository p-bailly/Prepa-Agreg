// Pyramide des chapitres (accueil) : Terminale à la base, agrégation au sommet, un trait par prérequis. Toucher un chapitre
// met en avant ce sur quoi il s'appuie et ce vers quoi il ouvre, et affiche sa fiche. Disposition calculée par
// src/contenu/graphe.ts ; teinte des nœuds selon l'avancement (une seule teinte, du clair au foncé, « à venir » en pointillés).
import { useLayoutEffect, useMemo, useRef, useState, type MouseEvent } from 'react'
import { GRAPHE } from '../contenu/chargement'
import { chapitresDe, disposerPyramide, liensChapitres } from '../contenu/graphe'
import { NOMS_NIVEAUX } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'
import { etatAffiche } from '../repetition/etat'

type Avancement = 'a-venir' | 'a-decouvrir' | 'en-cours' | 'acquis'

const NOMS_AVANCEMENT: Record<Avancement, string> = {
  'a-venir': 'À venir',
  'a-decouvrir': 'À découvrir',
  'en-cours': 'En cours',
  acquis: 'Acquis',
}

const ZOOMS = [1, 1.6, 2.4, 3.2]
const RAYON = 10

/** Nom court sur deux lignes au plus, coupé à l'espace (ou après le trait d'union) la plus proche du milieu. */
function lignes(nom: string): string[] {
  if (nom.length <= 11) return [nom]
  const coupures = [...nom.matchAll(/[ -]/g)].map((m) => m.index)
  if (!coupures.length) return [nom]
  const i = coupures.reduce((a, b) => (Math.abs(b - nom.length / 2) < Math.abs(a - nom.length / 2) ? b : a))
  return nom[i] === ' ' ? [nom.slice(0, i), nom.slice(i + 1)] : [nom.slice(0, i + 1), nom.slice(i + 1)]
}

export function Pyramide({ surChapitre }: { surChapitre: (id: string) => void }) {
  const { contenu, progression, reglages } = useDonnees()
  const liens = useMemo(() => liensChapitres(contenu, GRAPHE), [contenu])
  const p = useMemo(() => disposerPyramide(contenu, liens), [contenu, liens])
  const chapitres = useMemo(() => new Map(chapitresDe(contenu).map((ch) => [ch.id, ch])), [contenu])
  const themeDe = useMemo(() => new Map(contenu.themes.flatMap((t) => t.chapitres.map((ch) => [ch.id, t.titre]))), [contenu])
  const position = useMemo(() => new Map(p.noeuds.map((n) => [n.id, n])), [p])

  const avancement = useMemo(() => {
    const m = new Map<string, Avancement>()
    for (const id of chapitres.keys()) {
      const etats = contenu.cartes.filter((k) => k.chapitre === id).map((k) => etatAffiche(k, progression, reglages))
      m.set(
        id,
        !etats.length ? 'a-venir' : etats.every((e) => e === 'acquise') ? 'acquis' : etats.some((e) => e !== 'neuve') ? 'en-cours' : 'a-decouvrir',
      )
    }
    return m
  }, [chapitres, contenu, progression, reglages])

  const [choisi, setChoisi] = useState<string | null>(null)
  const [zoom, setZoom] = useState(0)
  const cadre = useRef<HTMLDivElement>(null)
  const svg = useRef<SVGSVGElement>(null)
  // Au changement de zoom, garder au centre le chapitre choisi (ou le même point qu'avant).
  const centreAvant = useRef<{ fx: number; fy: number } | null>(null)
  useLayoutEffect(() => {
    const c = cadre.current
    const f = centreAvant.current
    if (!c || !f) return
    c.scrollLeft = f.fx * c.scrollWidth - c.clientWidth / 2
    c.scrollTop = f.fy * c.scrollHeight - c.clientHeight / 2
    centreAvant.current = null
  }, [zoom])
  const changerZoom = (z: number) => {
    const c = cadre.current
    const n = choisi ? position.get(choisi) : undefined
    if (c)
      centreAvant.current = n
        ? { fx: n.x / p.largeur, fy: n.y / p.hauteur }
        : { fx: (c.scrollLeft + c.clientWidth / 2) / c.scrollWidth, fy: (c.scrollTop + c.clientHeight / 2) / c.scrollHeight }
    setZoom(Math.max(0, Math.min(ZOOMS.length - 1, z)))
  }

  // Un toucher choisit le chapitre le plus proche (les nœuds sont petits quand la pyramide tient dans l'écran).
  const toucher = (e: MouseEvent<SVGSVGElement>) => {
    const s = svg.current
    const m = s?.getScreenCTM()
    if (!s || !m) return
    const pt = new DOMPoint(e.clientX, e.clientY).matrixTransform(m.inverse())
    let meilleur: string | null = null
    let distance = 40
    for (const n of p.noeuds) {
      const d = Math.hypot(n.x - pt.x, n.y - pt.y)
      if (d < distance) [meilleur, distance] = [n.id, d]
    }
    setChoisi(meilleur === choisi ? null : meilleur)
  }

  const avant = choisi ? liens.filter((l) => l.vers === choisi).map((l) => l.de) : []
  const apres = choisi ? liens.filter((l) => l.de === choisi).map((l) => l.vers) : []
  const proches = new Set([choisi, ...avant, ...apres])
  const ch = choisi ? chapitres.get(choisi) : undefined
  const nbCartes = ch ? contenu.cartes.filter((k) => k.chapitre === ch.id).length : 0

  const courbe = (de: string, vers: string) => {
    const a = position.get(de)!
    const b = position.get(vers)!
    // Même rangée : un arc au-dessus, pour ne pas traverser les nœuds intermédiaires.
    if (Math.abs(a.y - b.y) < 1) {
      const h = Math.min(40, 12 + Math.abs(b.x - a.x) / 6)
      return `M${a.x},${a.y} C${a.x},${a.y - h} ${b.x},${b.y - h} ${b.x},${b.y}`
    }
    const my = (a.y + b.y) / 2
    return `M${a.x},${a.y} C${a.x},${my} ${b.x},${my} ${b.x},${b.y}`
  }

  return (
    <>
      <div className="pyramide-cadre" ref={cadre}>
        <svg
          ref={svg}
          className="pyramide"
          viewBox={`0 0 ${p.largeur} ${p.hauteur}`}
          style={{ width: `${ZOOMS[zoom] * 100}%` }}
          onClick={toucher}
          role="img"
          aria-label="Pyramide des chapitres, de la Terminale à l’agrégation, reliés par leurs prérequis"
        >
          {p.bandes.map((b, i) => (
            <g key={b.niveau}>
              <polygon
                className={i % 2 ? 'bande' : 'bande paire'}
                points={`${p.centre - b.demiHaut},${b.yHaut} ${p.centre + b.demiHaut},${b.yHaut} ${p.centre + b.demiBas},${b.yBas} ${p.centre - b.demiBas},${b.yBas}`}
              />
              <text className="niveau" x={8} y={(b.yHaut + b.yBas) / 2 + 4}>
                {NOMS_NIVEAUX[b.niveau]}
              </text>
            </g>
          ))}
          <g className={choisi ? 'liens estompes' : 'liens'}>
            {liens.map((l) => (
              <path key={`${l.de}>${l.vers}`} d={courbe(l.de, l.vers)} />
            ))}
          </g>
          {choisi && (
            <g className="liens-choisis">
              {avant.map((de) => (
                <path key={de} className="avant" d={courbe(de, choisi)} />
              ))}
              {apres.map((vers) => (
                <path key={vers} className="apres" d={courbe(choisi, vers)} />
              ))}
            </g>
          )}
          {p.noeuds.map((n) => {
            const nom = lignes(GRAPHE.courts[n.id] ?? n.id)
            return (
              <g
                key={n.id}
                className={`noeud-graphe ${avancement.get(n.id)}${choisi && !proches.has(n.id) ? ' loin' : ''}${n.id === choisi ? ' choisi' : ''}`}
                tabIndex={0}
                role="button"
                aria-label={`${chapitres.get(n.id)?.titre} (${NOMS_NIVEAUX[n.niveau]}, ${NOMS_AVANCEMENT[avancement.get(n.id)!]})`}
                onKeyDown={(e) => {
                  if (e.key === 'Enter' || e.key === ' ') {
                    e.preventDefault()
                    setChoisi(n.id === choisi ? null : n.id)
                  }
                }}
              >
                <title>{chapitres.get(n.id)?.titre}</title>
                {n.id === choisi && <circle className="halo" cx={n.x} cy={n.y} r={RAYON + 5} />}
                <circle cx={n.x} cy={n.y} r={RAYON} />
                <text x={n.x} y={n.y + RAYON + 12}>
                  {nom.map((l, i) => (
                    <tspan key={i} x={n.x} dy={i ? 11 : 0}>
                      {l}
                    </tspan>
                  ))}
                </text>
              </g>
            )
          })}
        </svg>
      </div>

      <div className="pyramide-outils">
        <ul className="legende-pyramide">
          {(Object.keys(NOMS_AVANCEMENT) as Avancement[]).map((a) => (
            <li key={a}>
              <svg width="14" height="14" viewBox="0 0 14 14" aria-hidden="true" className={`noeud-graphe ${a}`}>
                <circle cx="7" cy="7" r="5.5" />
              </svg>
              {NOMS_AVANCEMENT[a]}
            </li>
          ))}
        </ul>
        <div className="zoom">
          <button className="bouton compact" onClick={() => changerZoom(zoom - 1)} disabled={zoom === 0} aria-label="Dézoomer">
            −
          </button>
          <button className="bouton compact" onClick={() => changerZoom(zoom + 1)} disabled={zoom === ZOOMS.length - 1} aria-label="Zoomer">
            +
          </button>
        </div>
      </div>

      {ch ? (
        <div className="bloc fiche-chapitre">
          <span className="etiquette">
            {themeDe.get(ch.id)} · {NOMS_NIVEAUX[ch.niveau]} · {NOMS_AVANCEMENT[avancement.get(ch.id)!]}
          </span>
          <p className="titre-fiche">{ch.titre}</p>
          <p className="petit doux">
            {nbCartes ? `${nbCartes} carte${nbCartes > 1 ? 's' : ''} rédigée${nbCartes > 1 ? 's' : ''} sur ${ch.prevues}` : `${ch.prevues} notions prévues, pas encore rédigées`}
          </p>
          {avant.length > 0 && <LiensFiche titre="S’appuie sur" genre="avant" ids={avant} noms={chapitres} surChoix={setChoisi} />}
          {apres.length > 0 && <LiensFiche titre="Ouvre vers" genre="apres" ids={apres} noms={chapitres} surChoix={setChoisi} />}
          <span className="etiquette">Notions</span>
          <ul className="notions-fiche">
            {ch.notions.map((n) => (
              <li key={n}>{n}</li>
            ))}
          </ul>
          <button className="bouton" onClick={() => surChapitre(ch.id)}>
            Ouvrir le chapitre
          </button>
        </div>
      ) : (
        <p className="petit doux">Touchez un chapitre pour voir ce sur quoi il s’appuie et ce vers quoi il ouvre.</p>
      )}
    </>
  )
}

function LiensFiche({
  titre,
  genre,
  ids,
  noms,
  surChoix,
}: {
  titre: string
  genre: 'avant' | 'apres'
  ids: string[]
  noms: Map<string, { titre: string }>
  surChoix: (id: string) => void
}) {
  return (
    <>
      <span className="etiquette">
        <svg width="18" height="8" viewBox="0 0 18 8" aria-hidden="true" className={`trait-legende ${genre}`}>
          <line x1="1" y1="4" x2="17" y2="4" />
        </svg>
        {titre}
      </span>
      <div className="puces">
        {ids.map((id) => (
          <button key={id} className="puce" onClick={() => surChoix(id)}>
            {noms.get(id)?.titre ?? id}
          </button>
        ))}
      </div>
    </>
  )
}
