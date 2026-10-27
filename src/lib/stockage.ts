// Sauvegarde sur l'appareil (localStorage), protégée contre les navigateurs qui la refusent. Clés versionnées
// (`reglages-v1`, `progression-v1`) : changer un format = changer le suffixe ; la lecture fusionne avec les valeurs
// par défaut, pour qu'un champ ajouté plus tard ait sa valeur par défaut chez qui a déjà des données.

export function lire<T extends object>(cle: string, defaut: T): T {
  try {
    const t = localStorage.getItem(cle)
    if (t === null) return defaut
    const lu = JSON.parse(t) as unknown
    return lu && typeof lu === 'object' && !Array.isArray(lu) ? { ...defaut, ...(lu as Partial<T>) } : defaut
  } catch {
    return defaut
  }
}

export function ecrire(cle: string, valeur: unknown): void {
  try {
    localStorage.setItem(cle, JSON.stringify(valeur))
  } catch {
    /* stockage plein ou interdit : la valeur reste en mémoire pour cette session */
  }
}
