// Petits blocs communs du thème « Tonal ardoise » : encadrés typés, volets repliables, bandeau « Non garanti ».
import { useState, type ReactNode } from 'react'
import type { TypeCarte } from '../contenu/types'
import { Chevron } from './Icones'
import { Texte } from './Texte'

export type GenreEncadre = TypeCarte | 'exemple' | 'remarque'

const NOMS: Record<GenreEncadre, string> = {
  theoreme: 'Théorème',
  definition: 'Définition',
  methode: 'Méthode',
  exemple: 'Exemple',
  remarque: 'Remarque',
}

/**
 * Encadré d'un résultat : aplat de la couleur du type, nom du type en petit au-dessus du texte.
 * Ce nom porte le nom du type, suivi de `titre` ; `nom` le remplace entièrement.
 */
export function Encadre({ genre, titre, nom, children }: { genre: GenreEncadre; titre?: string; nom?: string; children: ReactNode }) {
  return (
    <section className={`encadre ${genre}`}>
      <span className="onglet">{nom ?? (titre ? `${NOMS[genre]} — ${titre}` : NOMS[genre])}</span>
      {children}
    </section>
  )
}

interface VoletProps {
  titre: string
  /** Petite mention à droite (« à revoir », « noté »). */
  marque?: string
  ouvert?: boolean
  surBascule?: (ouvert: boolean) => void
  children: ReactNode
}

/** Volet repliable : carte arrondie, chevron à droite (retourné quand il est ouvert). Contrôlé si `ouvert` est fourni. */
export function Volet({ titre, marque, ouvert, surBascule, children }: VoletProps) {
  const [interne, setInterne] = useState(false)
  const estOuvert = ouvert ?? interne
  const basculer = () => {
    if (ouvert === undefined) setInterne(!estOuvert)
    surBascule?.(!estOuvert)
  }
  return (
    <div className="volet">
      <button className="tete-volet" aria-expanded={estOuvert} onClick={basculer}>
        <span className="titre-volet">{titre}</span>
        {marque && <span className="marque">{marque}</span>}
        <span className="fleche">
          <Chevron />
        </span>
      </button>
      {estOuvert && <div className="corps-volet">{children}</div>}
    </div>
  )
}

/** Bandeau affiché quand une partie du contenu reste incertaine (champ `doute`). */
export function BandeauDoute({ doute }: { doute?: string }) {
  if (!doute) return null
  return (
    <div className="bandeau" role="note">
      <strong>Non garanti : </strong>
      <Texte texte={doute} enLigne />
    </div>
  )
}

export function Etiquette({ children }: { children: ReactNode }) {
  return <span className="etiquette">{children}</span>
}

/** Onglets internes d'un écran (Notions · Exercices · QCM…) : boutons collés en pilule, l'actif en aplat tonal. */
export function Segments<T extends string>({ options, valeur, surChange }: { options: { id: T; nom: string }[]; valeur: T; surChange: (v: T) => void }) {
  return (
    <div className="segments" role="tablist">
      {options.map((o) => (
        <button key={o.id} role="tab" aria-selected={o.id === valeur} onClick={() => surChange(o.id)}>
          {o.nom}
        </button>
      ))}
    </div>
  )
}
