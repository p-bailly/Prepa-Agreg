// Contrôle du contenu (content/) : identifiants, thèmes, niveaux, champs obligatoires, prérequis, trous, affirmations,
// exercices. Renvoie la liste des erreurs (vide si tout va bien). Utilisé par validation.test.ts ; voir CLAUDE.md,
// « Format du contenu ».

import { NIVEAUX, rangNiveau, type Carte, type Contenu, type Exercice } from './types'
import { trousAffiches, trousDe } from './texte'

const TYPES = ['theoreme', 'definition', 'methode']
const ID = /^[a-z0-9]+(-[a-z0-9]+)*$/

const texteNonVide = (x: unknown): x is string => typeof x === 'string' && x.trim() !== ''
const listeDeTextes = (x: unknown): x is string[] => Array.isArray(x) && x.every(texteNonVide)

export function validerContenu(c: Contenu): string[] {
  const err: string[] = []
  const themes = new Map(c.themes.map((t) => [t.id, t]))
  const chapitres = new Map<string, { theme: string; niveau: string }>()

  // Thèmes et chapitres.
  if (themes.size !== c.themes.length) err.push('themes.json : identifiant de thème en double')
  for (const t of c.themes) {
    if (!ID.test(t.id) || !texteNonVide(t.titre)) err.push(`thème ${t.id} : id ou titre invalide`)
    let rang = 0
    for (const ch of t.chapitres) {
      if (chapitres.has(ch.id)) err.push(`chapitre ${ch.id} en double`)
      if (!NIVEAUX.includes(ch.niveau)) err.push(`chapitre ${ch.id} : niveau inconnu ${ch.niveau}`)
      if (!listeDeTextes(ch.notions) || !ch.notions.length) err.push(`chapitre ${ch.id} : « notions » doit être une liste non vide`)
      else if (ch.prevues !== ch.notions.length) err.push(`chapitre ${ch.id} : « prevues » doit valoir le nombre de notions`)
      if (rangNiveau(ch.niveau) < rang) err.push(`thème ${t.id} : chapitres pas rangés par niveau croissant (${ch.id})`)
      rang = rangNiveau(ch.niveau)
      chapitres.set(ch.id, { theme: t.id, niveau: ch.niveau })
    }
  }

  // Cartes.
  const cartes = new Map<string, Carte>()
  for (const k of c.cartes) {
    const ou = `carte ${k.id}`
    if (!texteNonVide(k.id) || !ID.test(k.id)) err.push(`${ou} : identifiant invalide`)
    if (cartes.has(k.id)) err.push(`${ou} : identifiant en double`)
    cartes.set(k.id, k)
    if (!themes.has(k.theme)) err.push(`${ou} : thème inconnu ${k.theme}`)
    if (!NIVEAUX.includes(k.niveau)) err.push(`${ou} : niveau inconnu ${k.niveau}`)
    const ch = chapitres.get(k.chapitre)
    if (!ch) err.push(`${ou} : chapitre inconnu ${k.chapitre}`)
    else if (ch.theme !== k.theme || ch.niveau !== k.niveau) err.push(`${ou} : chapitre ${k.chapitre} d'un autre thème ou niveau`)
    if (!TYPES.includes(k.type)) err.push(`${ou} : type inconnu ${k.type}`)
    for (const champ of ['titre', 'recto', 'enonce'] as const) if (!texteNonVide(k[champ])) err.push(`${ou} : « ${champ} » manquant`)
    if (!Array.isArray(k.prerequis)) err.push(`${ou} : « prerequis » doit être une liste`)
    for (const champ of ['cours', 'demo', 'exemple', 'contre_exemple', 'doute'] as const)
      if (k[champ] !== undefined && !texteNonVide(k[champ])) err.push(`${ou} : « ${champ} » vide`)
    for (const champ of ['hypotheses', 'hypotheses_pieges'] as const)
      if (k[champ] !== undefined && !listeDeTextes(k[champ])) err.push(`${ou} : « ${champ} » doit être une liste de textes`)
    if (!!k.hypotheses?.length !== !!k.hypotheses_pieges?.length)
      err.push(`${ou} : « hypotheses » et « hypotheses_pieges » vont ensemble`)
    const toutes = [...(k.hypotheses ?? []), ...(k.hypotheses_pieges ?? [])]
    if (new Set(toutes).size !== toutes.length) err.push(`${ou} : hypothèses en double`)
    err.push(...validerTrous(k))
    if (k.affirmations !== undefined) {
      if (!Array.isArray(k.affirmations) || !k.affirmations.length) err.push(`${ou} : « affirmations » vide`)
      else {
        for (const a of k.affirmations) {
          if (!texteNonVide(a.texte) || typeof a.vrai !== 'boolean') err.push(`${ou} : affirmation mal formée`)
          if (a.vrai === false && !texteNonVide(a.contre_exemple)) err.push(`${ou} : affirmation fausse sans contre-exemple`)
        }
        if (!k.affirmations.some((a) => a.vrai === false)) err.push(`${ou} : le vrai/faux demande au moins une affirmation fausse`)
      }
    }
  }
  // Introductions de leçon : un chapitre existant, une seule par chapitre, un texte.
  const vues = new Set<string>()
  for (const l of c.lecons ?? []) {
    if (!chapitres.has(l.chapitre)) err.push(`leçon ${l.chapitre} : chapitre inconnu`)
    if (vues.has(l.chapitre)) err.push(`leçon ${l.chapitre} : en double`)
    vues.add(l.chapitre)
    if (!texteNonVide(l.introduction)) err.push(`leçon ${l.chapitre} : « introduction » manquante`)
  }

  // Prérequis : existants, jamais d'un niveau supérieur, sans cycle.
  for (const k of c.cartes) {
    for (const p of k.prerequis ?? []) {
      const q = cartes.get(p)
      if (!q) err.push(`carte ${k.id} : prérequis inconnu ${p}`)
      else if (rangNiveau(q.niveau) > rangNiveau(k.niveau)) err.push(`carte ${k.id} : prérequis ${p} d'un niveau supérieur`)
      if (p === k.id) err.push(`carte ${k.id} : prérequis d'elle-même`)
    }
  }
  const cycle = chercherCycle(cartes)
  if (cycle) err.push(`prérequis en cycle : ${cycle.join(' → ')}`)

  // Exercices.
  const ids = new Set<string>()
  const difficultes = new Map<string, Set<number>>()
  for (const e of c.exercices) {
    const ou = `exercice ${e.id}`
    if (!texteNonVide(e.id) || !ID.test(e.id)) err.push(`${ou} : identifiant invalide`)
    if (ids.has(e.id) || cartes.has(e.id)) err.push(`${ou} : identifiant en double`)
    ids.add(e.id)
    if (!cartes.has(e.carte)) err.push(`${ou} : carte inconnue ${e.carte}`)
    if (!(Number.isInteger(e.difficulte) && e.difficulte >= 1 && e.difficulte <= 5)) err.push(`${ou} : difficulté hors de 1..5`)
    for (const champ of ['titre', 'enonce', 'correction'] as const) if (!texteNonVide(e[champ])) err.push(`${ou} : « ${champ} » manquant`)
    if (!(typeof e.duree_min === 'number' && e.duree_min > 0)) err.push(`${ou} : « duree_min » invalide`)
    if (!listeDeTextes(e.indices)) err.push(`${ou} : « indices » doit être une liste de textes`)
    else if (e.difficulte >= 3 && e.indices.length < 2) err.push(`${ou} : au moins 2 indices à partir de la difficulté 3`)
    else if (e.indices.length > 3) err.push(`${ou} : 3 indices au plus`)
    err.push(...validerReponse(ou, e.reponse_finale))
    if (e.variantes !== undefined) {
      if (!Array.isArray(e.variantes) || !e.variantes.length) err.push(`${ou} : « variantes » vide`)
      else
        e.variantes.forEach((v, i) => {
          const ouV = `${ou}, variante ${i + 1}`
          for (const champ of ['enonce', 'correction'] as const) if (!texteNonVide(v[champ])) err.push(`${ouV} : « ${champ} » manquant`)
          if (v.enonce === e.enonce) err.push(`${ouV} : même énoncé que l'exercice`)
          if (v.titre !== undefined && !texteNonVide(v.titre)) err.push(`${ouV} : « titre » vide`)
          if (!!v.reponse_finale !== !!e.reponse_finale) err.push(`${ouV} : « reponse_finale » présente ici et pas dans l'exercice, ou l'inverse`)
          err.push(...validerReponse(ouV, v.reponse_finale))
          // Les indices d'un exercice citent souvent ses valeurs : chaque variante doit avoir les siens.
          if (e.indices?.length && v.indices === undefined) err.push(`${ouV} : « indices » propres à la variante manquants`)
          if (v.indices !== undefined && (!listeDeTextes(v.indices) || (e.difficulte >= 3 && v.indices.length < 2)))
            err.push(`${ouV} : « indices » invalides`)
        })
    }
    if (e.doute !== undefined && !texteNonVide(e.doute)) err.push(`${ou} : « doute » vide`)
    const d = difficultes.get(e.carte) ?? new Set<number>()
    d.add(e.difficulte)
    difficultes.set(e.carte, d)
  }
  // Questions de QCM écrites : rattachées à une notion du chapitre du fichier, bien formées, avec une explication.
  const idsQcm = new Set<string>()
  for (const fq of c.qcm ?? []) {
    if (!chapitres.has(fq.chapitre)) err.push(`qcm ${fq.chapitre} : chapitre inconnu`)
    if (!Array.isArray(fq.questions) || !fq.questions.length) err.push(`qcm ${fq.chapitre} : « questions » doit être une liste non vide`)
    for (const q of fq.questions ?? []) {
      const ou = `question ${q.id}`
      if (!texteNonVide(q.id) || !ID.test(q.id) || !q.id.startsWith(`${fq.chapitre}-q`)) err.push(`${ou} : identifiant attendu « ${fq.chapitre}-qN »`)
      if (idsQcm.has(q.id)) err.push(`${ou} : identifiant en double`)
      idsQcm.add(q.id)
      const k = cartes.get(q.carte)
      if (!k || k.chapitre !== fq.chapitre) err.push(`${ou} : carte ${q.carte} absente du chapitre ${fq.chapitre}`)
      if (!texteNonVide(q.enonce) || !texteNonVide(q.explication)) err.push(`${ou} : « enonce » et « explication » obligatoires`)
      if (q.genre === 'vrai-faux') {
        if (typeof q.vrai !== 'boolean' || q.propositions !== undefined || q.bonne !== undefined)
          err.push(`${ou} : un vrai-faux a « vrai » (booléen) et ni « propositions » ni « bonne »`)
      } else if (q.genre === 'choix') {
        if (!listeDeTextes(q.propositions) || q.propositions.length !== 4) err.push(`${ou} : 4 propositions`)
        else if (new Set(q.propositions.map((x) => x.trim())).size !== 4) err.push(`${ou} : propositions non distinctes`)
        if (!(Number.isInteger(q.bonne) && q.bonne! >= 0 && q.bonne! < 4)) err.push(`${ou} : « bonne » invalide`)
        if (q.vrai !== undefined) err.push(`${ou} : « vrai » réservé aux vrai-faux`)
      } else err.push(`${ou} : genre inconnu (vrai-faux ou choix)`)
      if (q.doute !== undefined && !texteNonVide(q.doute)) err.push(`${ou} : « doute » vide`)
    }
  }
  for (const [ch, n] of [...new Set((c.qcm ?? []).map((fq) => fq.chapitre))].map((ch) => [ch, (c.qcm ?? []).filter((fq) => fq.chapitre === ch).length] as const))
    if (n > 1) err.push(`qcm ${ch} : plusieurs fichiers`)

  // Par carte : de l'application directe (difficulté 1 ou 2) jusqu'à la question d'écrit (4 ou 5), 3 à 6 exercices.
  for (const k of c.cartes) {
    const d = [...(difficultes.get(k.id) ?? [])]
    const nb = c.exercices.filter((e) => e.carte === k.id).length
    if (nb < 3 || nb > 6) err.push(`carte ${k.id} : 3 à 6 exercices (${nb})`)
    if (d.length < 2) err.push(`carte ${k.id} : au moins 2 difficultés distinctes`)
    if (d.length && Math.min(...d) > 2) err.push(`carte ${k.id} : il manque un exercice de difficulté 1 ou 2`)
    if (d.length && Math.max(...d) < 4) err.push(`carte ${k.id} : il manque un exercice de difficulté 4 ou 5`)
  }
  return err
}

function validerTrous(k: Carte): string[] {
  const ou = `carte ${k.id}`
  if (k.enonce_trous === undefined && k.trous === undefined) return []
  if (!texteNonVide(k.enonce_trous) || !Array.isArray(k.trous) || !k.trous.length)
    return [`${ou} : « enonce_trous » et « trous » vont ensemble`]
  const err: string[] = []
  const vus = trousDe(k.enonce_trous)
  const attendus = k.trous.map((_, i) => i + 1)
  if (vus.length !== attendus.length || vus.some((n, i) => n !== attendus[i]))
    err.push(`${ou} : trous {{n}} mal numérotés (vus ${vus.join(',')}, attendus ${attendus.join(',')})`)
  // Un trou à l'intérieur d'une formule ne s'afficherait pas.
  if (trousAffiches(k.enonce_trous).length !== vus.length) err.push(`${ou} : un trou est à l'intérieur d'une formule`)
  k.trous.forEach((t, i) => {
    if (!texteNonVide(t.reponse) || !listeDeTextes(t.leurres)) err.push(`${ou} : trou ${i + 1} mal formé`)
    else if (t.leurres.length !== 3) err.push(`${ou} : trou ${i + 1} : exactement 3 leurres (${t.leurres.length})`)
    else if (new Set([t.reponse, ...t.leurres].map((x) => x.trim())).size !== 4)
      err.push(`${ou} : trou ${i + 1} : leurres non distincts de la réponse ou entre eux`)
  })
  return err
}

function validerReponse(ou: string, r: Exercice['reponse_finale']): string[] {
  if (r === undefined) return []
  if (!texteNonVide(r.question) || !listeDeTextes(r.propositions)) return [`${ou} : « reponse_finale » mal formée`]
  const err: string[] = []
  if (r.propositions.length !== 4) err.push(`${ou} : 4 propositions (${r.propositions.length})`)
  if (new Set(r.propositions.map((x) => x.trim())).size !== r.propositions.length) err.push(`${ou} : propositions non distinctes`)
  if (!(Number.isInteger(r.bonne) && r.bonne >= 0 && r.bonne < r.propositions.length)) err.push(`${ou} : « bonne » invalide`)
  return err
}

function chercherCycle(cartes: Map<string, Carte>): string[] | null {
  const etat = new Map<string, 1 | 2>()
  const pile: string[] = []
  const visiter = (id: string): string[] | null => {
    if (etat.get(id) === 2) return null
    if (etat.get(id) === 1) return [...pile.slice(pile.indexOf(id)), id]
    etat.set(id, 1)
    pile.push(id)
    for (const p of cartes.get(id)?.prerequis ?? []) {
      if (!cartes.has(p)) continue
      const c = visiter(p)
      if (c) return c
    }
    pile.pop()
    etat.set(id, 2)
    return null
  }
  for (const id of cartes.keys()) {
    const c = visiter(id)
    if (c) return c
  }
  return null
}

/** Tous les textes d'une carte ou d'un exercice qui peuvent contenir des formules, avec leur emplacement. */
export function textesDe(c: Contenu): { ou: string; texte: string }[] {
  const t: { ou: string; texte: string }[] = []
  const ajouter = (ou: string, x: unknown) => {
    if (typeof x === 'string') t.push({ ou, texte: x })
  }
  for (const k of c.cartes) {
    for (const champ of ['titre', 'recto', 'cours', 'enonce', 'enonce_trous', 'demo', 'exemple', 'contre_exemple', 'doute'] as const)
      ajouter(`${k.id}.${champ}`, k[champ])
    k.hypotheses?.forEach((h, i) => ajouter(`${k.id}.hypotheses[${i}]`, h))
    k.hypotheses_pieges?.forEach((h, i) => ajouter(`${k.id}.hypotheses_pieges[${i}]`, h))
    k.trous?.forEach((tr, i) => [tr.reponse, ...tr.leurres].forEach((x, j) => ajouter(`${k.id}.trous[${i}][${j}]`, x)))
    k.affirmations?.forEach((a, i) => {
      ajouter(`${k.id}.affirmations[${i}].texte`, a.texte)
      ajouter(`${k.id}.affirmations[${i}].contre_exemple`, a.contre_exemple)
      ajouter(`${k.id}.affirmations[${i}].justification`, a.justification)
    })
  }
  for (const fq of c.qcm ?? [])
    for (const q of fq.questions ?? []) {
      for (const champ of ['enonce', 'explication', 'doute'] as const) ajouter(`${q.id}.${champ}`, q[champ])
      q.propositions?.forEach((x, i) => ajouter(`${q.id}.propositions[${i}]`, x))
    }
  for (const l of c.lecons ?? []) {
    ajouter(`leçon ${l.chapitre}.introduction`, l.introduction)
    ajouter(`leçon ${l.chapitre}.doute`, l.doute)
  }
  for (const e of c.exercices) {
    for (const champ of ['titre', 'enonce', 'correction', 'doute'] as const) ajouter(`${e.id}.${champ}`, e[champ])
    e.indices?.forEach((x, i) => ajouter(`${e.id}.indices[${i}]`, x))
    const reponse = (ou: string, r: Exercice['reponse_finale']) => {
      if (!r) return
      ajouter(`${ou}.question`, r.question)
      r.propositions?.forEach((x, i) => ajouter(`${ou}.propositions[${i}]`, x))
    }
    reponse(e.id, e.reponse_finale)
    e.variantes?.forEach((v, k) => {
      const ou = `${e.id}.variantes[${k}]`
      ajouter(`${ou}.titre`, v.titre)
      ajouter(`${ou}.enonce`, v.enonce)
      ajouter(`${ou}.correction`, v.correction)
      v.indices?.forEach((x, i) => ajouter(`${ou}.indices[${i}]`, x))
      reponse(ou, v.reponse_finale)
    })
  }
  return t
}
