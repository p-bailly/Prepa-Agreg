// Navigation : quatre onglets (Accueil · Parcours · Stats · Réglages) et, par-dessus, une pile d'écrans (séance, carte,
// exercice, chapitre). Le bouton retour d'Android dépile.
import { useCallback, useEffect, useRef, useState } from 'react'
import { voletsDe } from './contenu/types'
import { FournisseurDonnees, useDonnees } from './etat/Donnees'
import { noterVolet, retenirMode } from './lib/progression'
import { natif } from './lib/notifications'
import { choisirMode } from './repetition/modes'
import { etatCarte } from './repetition/etat'
import { construireSeance } from './repetition/seance'
import { Accueil } from './pages/Accueil'
import { Chapitre } from './pages/Chapitre'
import { EcranCarte } from './pages/EcranCarte'
import { EcranExercice } from './pages/EcranExercice'
import { Parcours } from './pages/Parcours'
import { Reglages } from './pages/Reglages'
import { nouvelleSeance, Seance, type EtatSeance } from './pages/Seance'
import { Stats } from './pages/Stats'
import type { Mode } from './lib/reglages'

type Onglet = 'accueil' | 'parcours' | 'stats' | 'reglages'

type Ecran =
  | { k: 'seance' }
  | { k: 'chapitre'; id: string }
  | { k: 'carte'; id: string; mode: Mode; neuve: boolean }
  | { k: 'exercice'; id: string }

const ONGLETS: { id: Onglet; nom: string }[] = [
  { id: 'accueil', nom: 'Accueil' },
  { id: 'parcours', nom: 'Parcours' },
  { id: 'stats', nom: 'Stats' },
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
  const [seance, setSeance] = useState<EtatSeance | null>(null)
  const empiler = (e: Ecran) => {
    setPile((p) => [...p, e])
    window.scrollTo(0, 0)
  }
  const depiler = useCallback(() => {
    setPile((p) => p.slice(0, -1))
    window.scrollTo(0, 0)
  }, [])
  const changerSeance = useCallback((f: (s: EtatSeance) => EtatSeance) => setSeance((s) => (s ? f(s) : s)), [])

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

  const reviser = () => {
    setSeance(nouvelleSeance(construireSeance(d.contenu, d.progression, d.reglages, d.horloge(), d.aujourdhui()), d.horloge()))
    empiler({ k: 'seance' })
  }

  const ouvrirCarte = (id: string) => {
    const carte = d.cartesParId.get(id)
    if (!carte) return
    const neuve = etatCarte(carte, d.progression.volets) === 'neuve'
    const mode = choisirMode(carte, d.reglages.modesActifs, d.progression.dernierMode, neuve)
    d.changerProgression((p) => retenirMode(p, mode))
    empiler({ k: 'carte', id, mode, neuve })
  }

  const haut = pile[pile.length - 1]
  if (haut?.k === 'seance' && seance)
    return (
      <Seance
        seance={seance}
        changer={changerSeance}
        surExercice={(id) => empiler({ k: 'exercice', id })}
        surFin={() => {
          setSeance(null)
          depiler()
        }}
      />
    )
  if (haut?.k === 'exercice') {
    const e = d.exercicesParId.get(haut.id)
    if (e) return <EcranExercice key={e.id} exercice={e} surRetour={depiler} />
  }
  if (haut?.k === 'carte') {
    const carte = d.cartesParId.get(haut.id)
    if (carte)
      return (
        <EcranCarte
          key={carte.id}
          carte={carte}
          volets={voletsDe(carte)}
          neuve={haut.neuve}
          mode={haut.mode}
          surNote={(volet, note, principal) => {
            const t = d.horloge()
            d.changerProgression((p) => noterVolet(p, carte.id, volet, note, t, d.aujourdhui(), principal))
            if (principal) depiler()
          }}
          surExercice={(id) => empiler({ k: 'exercice', id })}
          surQuitter={depiler}
          libelleQuitter="Retour"
        />
      )
  }

  let contenu
  if (haut?.k === 'chapitre') contenu = <Chapitre id={haut.id} surCarte={ouvrirCarte} surRetour={depiler} />
  else if (onglet === 'accueil') contenu = <Accueil surReviser={reviser} />
  else if (onglet === 'parcours') contenu = <Parcours surChapitre={(id) => empiler({ k: 'chapitre', id })} />
  else if (onglet === 'stats') contenu = <Stats />
  else contenu = <Reglages />

  return (
    <div className="appli avec-onglets" style={onglet === 'accueil' && !haut ? { ['--pied' as string]: '100px' } : undefined}>
      {contenu}
      <nav className="onglets" aria-label="Onglets">
        {ONGLETS.map((o) => (
          <button
            key={o.id}
            aria-current={o.id === onglet && !haut ? 'page' : undefined}
            onClick={() => {
              setPile([])
              setOnglet(o.id)
              window.scrollTo(0, 0)
            }}
          >
            {o.nom}
          </button>
        ))}
      </nav>
    </div>
  )
}
