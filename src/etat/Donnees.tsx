// Données de l'appli partagées par les écrans : contenu embarqué, réglages et progression persistés, horloge.
import { createContext, useCallback, useContext, useEffect, useMemo, useState, type ReactNode } from 'react'
import { CONTENU } from '../contenu/chargement'
import type { Carte, Contenu, Exercice } from '../contenu/types'
import { ecrireProgression, lireProgression, type Progression } from '../lib/progression'
import { ecrireReglages, estSombre, lireReglages, type Reglages } from '../lib/reglages'
import { horlogeSysteme, jourLocal, type Horloge } from '../lib/temps'

interface Donnees {
  contenu: Contenu
  cartesParId: Map<string, Carte>
  exercicesParCarte: Map<string, Exercice[]>
  exercicesParId: Map<string, Exercice>
  reglages: Reglages
  changerReglages: (f: (r: Reglages) => Reglages) => void
  progression: Progression
  changerProgression: (f: (p: Progression) => Progression) => void
  /** Remplace toute la progression (import d'une sauvegarde). */
  remplacer: (p: Progression, r?: Reglages) => void
  horloge: Horloge
  /** Jour civil local de l'instant présent. */
  aujourdhui: () => string
}

const Contexte = createContext<Donnees | null>(null)

export function FournisseurDonnees({ children, horloge = horlogeSysteme }: { children: ReactNode; horloge?: Horloge }) {
  const [reglages, setReglages] = useState(lireReglages)
  const [progression, setProgression] = useState(lireProgression)

  useEffect(() => ecrireReglages(reglages), [reglages])
  useEffect(() => ecrireProgression(progression), [progression])
  useEffect(() => {
    document.documentElement.dataset.theme = estSombre(reglages) ? 'dark' : 'light'
  }, [reglages])

  const changerReglages = useCallback((f: (r: Reglages) => Reglages) => setReglages(f), [])
  const changerProgression = useCallback((f: (p: Progression) => Progression) => setProgression(f), [])
  const remplacer = useCallback((p: Progression, r?: Reglages) => {
    setProgression(p)
    if (r) setReglages(r)
  }, [])

  const index = useMemo(() => {
    const cartesParId = new Map(CONTENU.cartes.map((c) => [c.id, c]))
    const exercicesParId = new Map(CONTENU.exercices.map((e) => [e.id, e]))
    const exercicesParCarte = new Map<string, Exercice[]>()
    for (const e of CONTENU.exercices) exercicesParCarte.set(e.carte, [...(exercicesParCarte.get(e.carte) ?? []), e])
    for (const l of exercicesParCarte.values()) l.sort((a, b) => a.difficulte - b.difficulte)
    return { cartesParId, exercicesParId, exercicesParCarte }
  }, [])

  const valeur = useMemo<Donnees>(
    () => ({
      contenu: CONTENU,
      ...index,
      reglages,
      changerReglages,
      progression,
      changerProgression,
      remplacer,
      horloge,
      aujourdhui: () => jourLocal(horloge()),
    }),
    [index, reglages, changerReglages, progression, changerProgression, remplacer, horloge],
  )
  return <Contexte.Provider value={valeur}>{children}</Contexte.Provider>
}

export function useDonnees(): Donnees {
  const d = useContext(Contexte)
  if (!d) throw new Error('useDonnees hors de FournisseurDonnees')
  return d
}
