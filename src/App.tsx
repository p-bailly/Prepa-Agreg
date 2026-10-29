// Navigation : trois onglets (Accueil · Parcours · Réglages) et, par-dessus, une pile d'écrans (chapitre, notion, QCM,
// carte, exercice). Tous les écrans de la pile restent montés (seul le dernier est visible) : au retour, chacun est tel
// qu'on l'a laissé. Le bouton retour d'Android dépile.
import { useCallback, useEffect, useRef, useState, type ReactNode } from 'react'
import { cleQcm, type CibleQcm } from './chapitres/qcm'
import { voletsDe, type Niveau } from './contenu/types'
import { FournisseurDonnees, useDonnees } from './etat/Donnees'
import { noterVolet, retenirMode } from './lib/progression'
import { natif } from './lib/notifications'
import { choisirMode } from './repetition/modes'
import { etatCarte } from './repetition/etat'
import { IconeOnglet } from './composants/Icones'
import { Accueil } from './pages/Accueil'
import { Chapitre } from './pages/Chapitre'
import { EcranCarte } from './pages/EcranCarte'
import { EcranExercice } from './pages/EcranExercice'
import { EcranNotion } from './pages/EcranNotion'
import { EcranQcm } from './pages/EcranQcm'
import { Parcours } from './pages/Parcours'
import { Reglages } from './pages/Reglages'
import type { Mode } from './lib/reglages'

type Onglet = 'accueil' | 'parcours' | 'reglages'

type Ecran =
  | { k: 'chapitre'; id: string }
  | { k: 'notion'; id: string }
  | { k: 'qcm'; cible: CibleQcm }
  | { k: 'carte'; id: string; mode: Mode; neuve: boolean }
  | { k: 'exercice'; id: string }

/** Écrans affichés avec la barre d'onglets ; les autres prennent tout l'écran. */
const AVEC_ONGLETS: Ecran['k'][] = ['chapitre', 'notion']

const ONGLETS: { id: Onglet; nom: string }[] = [
  { id: 'accueil', nom: 'Accueil' },
  { id: 'parcours', nom: 'Parcours' },
  { id: 'reglages', nom: 'Réglages' },
]

export default function App() {
  return (
    <FournisseurDonnees>
      <Navigation />
    </FournisseurDonnees>
  )
}

function Navigation() {
  const d = useDonnees()
  const [onglet, setOnglet] = useState<Onglet>('accueil')
  const [pile, setPile] = useState<Ecran[]>([])
  // Position de défilement de chaque écran de la pile, rendue au retour.
  const defilements = useRef<number[]>([])
  const empiler = (e: Ecran) => {
    defilements.current.push(window.scrollY)
    setPile((p) => [...p, e])
    window.scrollTo(0, 0)
  }
  const depiler = useCallback(() => {
    const y = defilements.current.pop() ?? 0
    setPile((p) => p.slice(0, -1))
    requestAnimationFrame(() => window.scrollTo(0, y))
  }, [])

  // Bouton retour d'Android : dépile, sinon revient à l'accueil, sinon laisse l'appli se fermer.
  const etat = useRef({ pile, onglet })
  etat.current = { pile, onglet }
  useEffect(() => {
    if (!natif) return
    let retrait: (() => void) | undefined
    void import('@capacitor/app').then(({ App: AppNatif }) =>
      AppNatif.addListener('backButton', () => {
        const { pile: p, onglet: o } = etat.current
        if (p.length) depiler()
        else if (o !== 'accueil') setOnglet('accueil')
        else void AppNatif.exitApp()
      }).then((h) => (retrait = () => void h.remove())),
    )
    return () => retrait?.()
  }, [depiler])

  // Même chose avec le bouton retour du navigateur (version web, GitHub Pages) : une entrée d'historique « sentinelle »
  // est posée tant qu'on n'est pas sur l'accueil sans écran empilé ; la retirer (retour) fait reculer l'appli d'un cran.
  const sentinelle = useRef(false)
  useEffect(() => {
    if (natif) return
    const surRetour = () => {
      sentinelle.current = false
      const { pile: p, onglet: o } = etat.current
      if (p.length) depiler()
      else if (o !== 'accueil') setOnglet('accueil')
    }
    window.addEventListener('popstate', surRetour)
    return () => window.removeEventListener('popstate', surRetour)
  }, [depiler])
  useEffect(() => {
    if (natif || sentinelle.current || (!pile.length && onglet === 'accueil')) return
    history.pushState(null, '')
    sentinelle.current = true
  }, [pile, onglet])

  const ouvrirCarte = (id: string) => {
    const carte = d.cartesParId.get(id)
    if (!carte) return
    const neuve = etatCarte(carte, d.progression) === 'neuve'
    const mode = choisirMode(carte, d.reglages.modesActifs, d.progression.dernierMode, neuve)
    d.changerProgression((p) => retenirMode(p, mode))
    empiler({ k: 'carte', id, mode, neuve })
  }
  const surChapitre = (id: string) => empiler({ k: 'chapitre', id })
  const surExercice = (id: string) => empiler({ k: 'exercice', id })

  const rendre = (e: Ecran): ReactNode => {
    if (e.k === 'chapitre')
      return (
        <Chapitre
          id={e.id}
          surNotion={(id) => empiler({ k: 'notion', id })}
          surExercice={surExercice}
          surQcm={() => empiler({ k: 'qcm', cible: { chapitre: e.id } })}
          surRetour={depiler}
        />
      )
    if (e.k === 'notion') {
      const carte = d.cartesParId.get(e.id)
      return carte && <EcranNotion carte={carte} surExercice={surExercice} surMemoriser={() => ouvrirCarte(carte.id)} surRetour={depiler} />
    }
    if (e.k === 'qcm') return <EcranQcm cible={e.cible} surFin={depiler} />
    if (e.k === 'exercice') {
      const exercice = d.exercicesParId.get(e.id)
      return exercice && <EcranExercice exercice={exercice} surRetour={depiler} />
    }
    const carte = d.cartesParId.get(e.id)
    return (
      carte && (
        <EcranCarte
          carte={carte}
          volets={voletsDe(carte)}
          neuve={e.neuve}
          mode={e.mode}
          surNote={(volet, note, principal) => {
            const t = d.horloge()
            d.changerProgression((p) => noterVolet(p, carte.id, volet, note, t, d.aujourdhui(), principal))
            if (principal) depiler()
          }}
          surExercice={surExercice}
          surQuitter={depiler}
          libelleQuitter="Retour"
        />
      )
    )
  }

  const haut = pile[pile.length - 1]
  const onglets = !haut || AVEC_ONGLETS.includes(haut.k)
  let base
  if (onglet === 'accueil')
    base = (
      <Accueil
        surChapitre={surChapitre}
        surQcmNiveau={(niveau: Niveau, sorte: 'qcm' | 'vrai-faux') => empiler({ k: 'qcm', cible: { niveau, sorte } })}
      />
    )
  else if (onglet === 'parcours') base = <Parcours surChapitre={surChapitre} />
  else base = <Reglages />

  return (
    <>
      <div className="appli avec-onglets" style={haut ? { display: 'none' } : undefined}>
        {base}
      </div>
      {pile.map((e, i) => (
        <div
          key={`${i}-${e.k}-${'id' in e ? e.id : cleQcm(e.cible)}`}
          className={AVEC_ONGLETS.includes(e.k) ? 'appli avec-onglets' : undefined}
          style={i < pile.length - 1 ? { display: 'none' } : undefined}
        >
          {rendre(e)}
        </div>
      ))}
      {onglets && (
        <nav className="onglets" aria-label="Onglets">
          {ONGLETS.map((o) => (
            <button
              key={o.id}
              aria-current={o.id === onglet && !haut ? 'page' : undefined}
              onClick={() => {
                defilements.current = []
                setPile([])
                setOnglet(o.id)
                window.scrollTo(0, 0)
              }}
            >
              <span className="indicateur">
                <IconeOnglet onglet={o.id} />
              </span>
              {o.nom}
            </button>
          ))}
        </nav>
      )}
    </>
  )
}
