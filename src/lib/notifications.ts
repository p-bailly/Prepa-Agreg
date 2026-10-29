// Plateforme (appli Android ou navigateur) et annulation de l'ancien rappel quotidien, retiré le 2026-10-29 : une
// version précédente a pu le programmer auprès d'Android, où il partirait même appli fermée.
import { Capacitor } from '@capacitor/core'
import { LocalNotifications } from '@capacitor/local-notifications'

export const natif = Capacitor.isNativePlatform()
const ID_RAPPEL = 1

export async function annulerRappel(): Promise<void> {
  if (!natif) return
  try {
    await LocalNotifications.cancel({ notifications: [{ id: ID_RAPPEL }] })
  } catch {
    /* rien à annuler */
  }
}
