// Chapitre : en tête son numéro, son thème et son niveau, puis trois onglets. Notions : introduction du chapitre et
// notions numérotées (chacune ouvre son cours et ses exercices). Exercices : tous les exercices du chapitre, les
// exercices (difficulté 1 à 3) puis les problèmes (4 et 5). QCM : questions tirées du chapitre, meilleur résultat.
import { useMemo, useState } from 'react'
import { avancementNotion, chapitresNumerotes, meilleurQcm, notionsDe } from '../chapitres/plan'
import { questionsChapitre } from '../chapitres/qcm'
import { BandeauDoute, Segments, Volet } from '../composants/Blocs'
import { IconeTheme } from '../composants/Icones'
import { Ligne, LigneExercice, NOMS_AVANCEMENT } from '../composants/Lignes'
import { Texte } from '../composants/Texte'
import { NOMS_NIVEAUX, type Exercice } from '../contenu/types'
import { useDonnees } from '../etat/Donnees'
import { etatExercice } from '../exercices/niveaux'

type Onglet = 'notions' | 'exercices' | 'qcm'

const NOMS_TYPES = { theoreme: 'Théorème', definition: 'Définition', methode: 'Méthode' } as const

export function Chapitre({
  id,
  surNotion,
  surExercice,
  surQcm,
  surRetour,
}: {
  id: string
  surNotion: (id: string) => void
  surExercice: (id: string) => void
  surQcm: () => void
  surRetour: () => void
}) {
  const { contenu, progression, exercicesParCarte } = useDonnees()
  const [onglet, setOnglet] = useState<Onglet>('notions')
  const entree = useMemo(() => chapitresNumerotes(contenu).find((x) => x.chapitre.id === id), [contenu, id])
  const notions = useMemo(() => notionsDe(contenu, id), [contenu, id])
  const nbQuestions = useMemo(() => questionsChapitre(contenu, id).length, [contenu, id])
  if (!entree) return null
  const { numero, chapitre, theme } = entree
  const intro = contenu.lecons?.find((l) => l.chapitre === id)
  const numeroNotion = new Map(notions.map((k, i) => [k.id, `${numero}.${i + 1}`]))
  const exercices = notions.flatMap((k) => exercicesParCarte.get(k.id) ?? [])
  const parDifficulte = (l: Exercice[]) => [...l].sort((a, b) => a.difficulte - b.difficulte)
  const qcm = meilleurQcm(progression, id)
  const faits = progression.qcm[id] ?? []

  const listeExercices = (l: Exercice[]) => (
    <ul className="liste">
      {parDifficulte(l).map((e) => (
        <li key={e.id}>
          <LigneExercice exercice={e} etat={etatExercice(progression.exercices[e.id])} avant={`Notion ${numeroNotion.get(e.carte)}`} surExercice={surExercice} />
        </li>
      ))}
    </ul>
  )

  return (
    <main className={`ecran th-${theme.id}`}>
      <button className="retour" onClick={surRetour}>
        ← Retour
      </button>
      <span className="etiquette-theme">
        <span className="rond">
          <IconeTheme theme={theme.id} />
        </span>
        Chapitre {numero} · {theme.titre} · {NOMS_NIVEAUX[chapitre.niveau]}
      </span>
      <h1 className="titre-ecran">{chapitre.titre}</h1>
      <Segments
        options={[
          { id: 'notions', nom: 'Notions' },
          { id: 'exercices', nom: 'Exercices' },
          { id: 'qcm', nom: 'QCM' },
        ]}
        valeur={onglet}
        surChange={setOnglet}
      />

      {onglet === 'notions' && (
        <>
          {intro && (
            <Volet titre="Introduction du chapitre">
              <div style={{ paddingTop: 12 }}>
                <BandeauDoute doute={intro.doute} />
                <Texte texte={intro.introduction} />
              </div>
            </Volet>
          )}
          <ul className="liste bloc-liste">
            {notions.map((k) => {
              const ex = exercicesParCarte.get(k.id) ?? []
              const a = avancementNotion(k, ex, progression)
              return (
                <li key={k.id}>
                  <Ligne
                    avancement={a}
                    numero={numeroNotion.get(k.id)}
                    titre={k.titre}
                    dessous={`${NOMS_TYPES[k.type]} · cours et ${ex.length} exercice${ex.length > 1 ? 's' : ''} · ${NOMS_AVANCEMENT[a]}`}
                    surClic={() => surNotion(k.id)}
                  />
                </li>
              )
            })}
          </ul>
          {chapitre.notions.length > notions.length && (
            <p className="petit doux">
              {notions.length} notions rédigées sur {chapitre.notions.length} prévues.
            </p>
          )}
        </>
      )}

      {onglet === 'exercices' && (
        <>
          <h2 className="sous-titre">Exercices ({exercices.filter((e) => e.difficulte <= 3).length})</h2>
          <p className="petit doux">Difficulté 1 à 3 : application directe, calcul standard, classique en plusieurs étapes.</p>
          <div className="bloc-liste">{listeExercices(exercices.filter((e) => e.difficulte <= 3))}</div>
          <h2 className="sous-titre">Problèmes ({exercices.filter((e) => e.difficulte >= 4).length})</h2>
          <p className="petit doux">Difficulté 4 et 5 : questions d’écrit enchaînées, idées moins immédiates.</p>
          {exercices.some((e) => e.difficulte >= 4) ? (
            <div className="bloc-liste">{listeExercices(exercices.filter((e) => e.difficulte >= 4))}</div>
          ) : (
            <p className="doux">Pas encore de problème pour ce chapitre.</p>
          )}
        </>
      )}

      {onglet === 'qcm' && (
        <div className="bloc">
          <p>
            <strong>{Math.min(10, nbQuestions)} questions</strong> tirées au hasard parmi {nbQuestions} : vrai ou faux, énoncés à
            compléter et questions d’exercices, avec la justification après chaque réponse.
          </p>
          {qcm ? (
            <p className="petit doux">
              Meilleur résultat : {qcm.bonnes} / {qcm.total} · fait {faits.length} fois.
            </p>
          ) : (
            <p className="petit doux">Pas encore fait.</p>
          )}
          <button className="bouton-principal" disabled={!nbQuestions} onClick={surQcm}>
            Lancer le QCM
          </button>
        </div>
      )}
    </main>
  )
}
