// Écran d'un exercice : énoncé, recherche sur papier, indices à débloquer un par un, réponse à choisir parmi quatre
// propositions (exercice à résultat), correction, puis auto-évaluation (Réussi, En partie, Raté). L'historique est gardé
// par exercice ; les exercices ne sont pas replanifiés.
import { useMemo, useState, type CSSProperties } from 'react'
import { BandeauDoute, Volet } from '../composants/Blocs'
import { Signaler } from '../composants/Signalement'
import { Texte } from '../composants/Texte'
import { nombreVariantes, NOMS_NIVEAUX, variante, type Exercice } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'
import { etatExercice, varianteSuivante } from '../exercices/niveaux'
import { melanger } from '../lib/hasard'
import { enregistrerTentative, type ResultatExercice } from '../lib/progression'

const basculer = (s: Set<number>, i: number, dans: boolean) => {
  const n = new Set(s)
  if (dans) n.add(i)
  else n.delete(i)
  return n
}

const RESULTATS: { r: ResultatExercice; nom: string }[] = [
  { r: 'reussi', nom: 'Réussi' },
  { r: 'partiel', nom: 'En partie' },
  { r: 'rate', nom: 'Raté' },
]

export function EcranExercice({ exercice: base, surRetour }: { exercice: Exercice; surRetour: () => void }) {
  const { progression } = useDonnees()
  const nombre = nombreVariantes(base)
  // Variante choisie à l'ouverture (jamais faite, sinon la plus ancienne) ; « Autres valeurs » passe à la suivante.
  const [k, setK] = useState(() => varianteSuivante(nombre, progression.exercices[base.id]))
  return <Deroule key={k} exercice={variante(base, k)} k={k} nombre={nombre} autre={() => setK((k + 1) % nombre)} surRetour={surRetour} />
}

function Deroule({
  exercice,
  k,
  nombre,
  autre,
  surRetour,
}: {
  exercice: Exercice
  k: number
  nombre: number
  autre: () => void
  surRetour: () => void
}) {
  const { cartesParId, contenu, progression, changerProgression, horloge, aujourdhui } = useDonnees()
  const carte = cartesParId.get(exercice.carte)
  const theme = contenu.themes.find((t) => t.id === carte?.theme)
  const historique = progression.exercices[exercice.id] ?? []
  const [indices, setIndices] = useState(0)
  const [replies, setReplies] = useState<Set<number>>(new Set())
  const [correction, setCorrection] = useState(false)
  const [choix, setChoix] = useState<number | null>(null)
  const [choisir, setChoisir] = useState(false)
  const rf = exercice.reponse_finale
  const propositions = useMemo(
    () => (rf ? melanger(rf.propositions.map((t, i) => ({ t, bonne: i === rf.bonne }))) : []),
    [rf],
  )

  const evaluer = (resultat: ResultatExercice) => {
    const t = horloge()
    changerProgression((p) =>
      enregistrerTentative(p, exercice.id, {
        date: t,
        jour: aujourdhui(),
        resultat,
        indices,
        variante: k,
        ...(choix !== null ? { reponse: propositions[choix].bonne ? 'juste' : 'fausse' } : {}),
      }),
    )
    surRetour()
  }

  const style = { '--pied': correction ? '110px' : '100px' } as CSSProperties
  const derniere = historique[historique.length - 1]

  return (
    <div className="appli" style={style}>
      <main className="ecran">
        <div className="barre-haut">
          <button className="retour" onClick={surRetour}>
            ← Retour à la carte
          </button>
          <span className="avancement">Difficulté {exercice.difficulte} / 5</span>
        </div>
        <header className="tete-carte">
          <span className="sur-titre">
            Exercice · {theme?.titre ?? carte?.theme} · {carte ? NOMS_NIVEAUX[carte.niveau] : ''}
          </span>
          <h1 className="titre-carte">
            <Texte texte={exercice.titre} enLigne />
          </h1>
          <p className="petit">Notion : {carte?.titre}</p>
        </header>
        <BandeauDoute doute={exercice.doute} />
        <div className="bloc">
          <Texte texte={exercice.enonce} />
        </div>
        {nombre > 1 && (
          <p className="petit doux">
            Variante {k + 1} / {nombre} (mêmes questions, autres valeurs).{' '}
            <button className="lien" style={{ padding: 0, minHeight: 0 }} onClick={autre}>
              Autres valeurs
            </button>
          </p>
        )}
        <p className="consigne">
          Cherchez sur papier (environ {exercice.duree_min} min).
          {derniere && ` Déjà fait ${historique.length} fois ; ${etatExercice(historique) === 'reussi' ? 'réussi' : 'pas encore réussi'}.`}
        </p>

        {exercice.indices.slice(0, indices).map((t, i) => (
          <Volet key={i} titre={`Indice ${i + 1}`} ouvert={!replies.has(i)} surBascule={(o) => setReplies((s) => basculer(s, i, !o))}>
            <div style={{ paddingTop: 12 }}>
              <Texte texte={t} />
            </div>
          </Volet>
        ))}
        {indices < exercice.indices.length && !correction && (
          <button className="bouton" onClick={() => setIndices(indices + 1)}>
            Indice {indices + 1} / {exercice.indices.length}
          </button>
        )}

        {rf && (
          <div style={{ marginTop: 10 }}>
            {!choisir ? (
              <button className="bouton" onClick={() => setChoisir(true)}>
                Choisir la réponse
              </button>
            ) : (
              <>
                <span className="etiquette">Réponse</span>
                <p>
                  <Texte texte={rf.question} enLigne />
                </p>
                <ul className="liste-choix">
                  {propositions.map((p, i) => {
                    const pris = choix === i
                    const classe = choix === null ? '' : p.bonne ? ' juste' : pris ? ' fausse' : ''
                    return (
                      <li key={i}>
                        <button className={`case${classe}`} aria-pressed={pris} disabled={choix !== null} onClick={() => setChoix(i)}>
                          <span className="coche" aria-hidden="true">
                            {pris ? '●' : ''}
                          </span>
                          <span className="contenu">
                            <Texte texte={p.t} enLigne />
                            {choix !== null && p.bonne && <span className="verdict juste">✓ Bonne réponse</span>}
                            {pris && !p.bonne && <span className="verdict fausse">✗ Mauvaise réponse</span>}
                          </span>
                        </button>
                      </li>
                    )
                  })}
                </ul>
              </>
            )}
          </div>
        )}

        <div style={{ marginTop: 12 }}>
          <Volet titre="Correction" ouvert={correction} surBascule={setCorrection}>
            <div style={{ paddingTop: 12 }}>
              <Texte texte={exercice.correction} />
            </div>
          </Volet>
        </div>
        <Signaler cible={exercice.id} genre="exercice" ou={correction ? 'correction' : undefined} />
      </main>

      <div className="pied">
        <div className="pied-interieur">
          {!correction ? (
            <button className="bouton-principal" onClick={() => setCorrection(true)}>
              Voir la correction
            </button>
          ) : (
            <div className="rangee-boutons trois">
              {RESULTATS.map((x) => (
                <button key={x.r} className="bouton" onClick={() => evaluer(x.r)}>
                  {x.nom}
                </button>
              ))}
            </div>
          )}
        </div>
      </div>
    </div>
  )
}
