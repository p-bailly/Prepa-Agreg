// Séance courte : enchaîne les cartes de la séance du jour (au plus N cartes ou N minutes), alterne les modes, remet
// en fin de séance une carte notée « À revoir », puis affiche le bilan.
import { useEffect } from 'react'
import type { Volet } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'
import { noterVolet, retenirMode } from '../lib/progression'
import { choisirMode } from '../repetition/modes'
import { NOMS_NOTES, NOTES, type Note } from '../repetition/planification'
import { avancer, position, type EtatSeance } from '../repetition/deroulement'
import { EcranCarte } from './EcranCarte'

export function Seance({
  seance,
  changer,
  surExercice,
  surFin,
}: {
  seance: EtatSeance
  changer: (f: (s: EtatSeance) => EtatSeance) => void
  surExercice: (id: string) => void
  surFin: () => void
}) {
  const { cartesParId, progression, changerProgression, reglages, horloge, aujourdhui } = useDonnees()
  const element = seance.elements[seance.index]
  const carte = element && cartesParId.get(element.carte)

  // Mode de la carte, choisi une fois à son arrivée (puis gardé dans l'état de séance) et retenu pour l'alternance.
  const modePropose =
    !seance.fin && carte && element ? choisirMode(carte, reglages.modesActifs, progression.dernierMode, element.neuve) : null
  const mode = seance.mode ?? modePropose
  useEffect(() => {
    if (seance.mode !== null || modePropose === null) return
    changer((s) => (s.mode === null ? { ...s, mode: modePropose } : s))
    changerProgression((p) => retenirMode(p, modePropose))
  }, [seance.mode, modePropose, changer, changerProgression])

  if (seance.fin || !carte || !element) return <Bilan seance={seance} surFin={surFin} />

  const noter = (volet: Volet, note: Note, principal: boolean) => {
    const t = horloge()
    const premiere = !seance.vues.includes(carte.id)
    changerProgression((p) => noterVolet(p, carte.id, volet, note, t, aujourdhui(), principal && premiere))
    if (!principal) return
    changer((s) => avancer(s, volet, note, t, reglages.seanceMinutes))
  }

  return (
    <EcranCarte
      key={`${seance.index}-${carte.id}`}
      carte={carte}
      position={position(seance)}
      volets={element.volets}
      neuve={element.neuve}
      mode={mode ?? 'tete'}
      surNote={noter}
      surExercice={surExercice}
      surQuitter={() => changer((s) => ({ ...s, fin: true }))}
      libelleQuitter="Terminer la séance"
    />
  )
}

function Bilan({ seance, surFin }: { seance: EtatSeance; surFin: () => void }) {
  const { horloge } = useDonnees()
  const minutes = Math.max(1, Math.round((horloge() - seance.debut) / 60_000))
  const vide = seance.elements.length === 0
  return (
    <div className="appli" style={{ ['--pied' as string]: '100px' }}>
      <main className="ecran">
        <h1 className="titre-ecran">Bilan</h1>
        {vide ? (
          <div className="bloc">
            <p>Rien à réviser pour l’instant : pas de carte due, et le plafond de cartes neuves du jour est atteint (ou les filtres ne laissent rien passer).</p>
          </div>
        ) : (
          <>
            <div className="trois-colonnes">
              <div>
                <div className="grand-nombre">{seance.vues.length}</div>
                <div className="petit">cartes</div>
              </div>
              <div>
                <div className="grand-nombre">{minutes}</div>
                <div className="petit">minutes</div>
              </div>
              <div>
                <div className="grand-nombre">{seance.notes.filter((x) => x >= 2).length}</div>
                <div className="petit">bien ou facile</div>
              </div>
            </div>
            <span className="etiquette">Répartition des notes</span>
            <div className="bloc">
              {NOTES.map((n) => {
                const k = seance.notes.filter((x) => x === n).length
                return (
                  <div key={n} className="jauge">
                    <div className="ligne">
                      <span>{NOMS_NOTES[n]}</span>
                      <span>{k}</span>
                    </div>
                    <div className="barre">
                      <span style={{ width: `${seance.notes.length ? (100 * k) / seance.notes.length : 0}%` }} />
                    </div>
                  </div>
                )
              })}
            </div>
            {seance.index < seance.elements.length && (
              <p className="doux">Séance arrêtée avant la fin : les cartes restantes attendent la prochaine séance.</p>
            )}
          </>
        )}
      </main>
      <div className="pied">
        <div className="pied-interieur">
          <button className="bouton-principal" onClick={surFin}>
            Terminer
          </button>
        </div>
      </div>
    </div>
  )
}
