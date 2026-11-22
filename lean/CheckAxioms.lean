import DFL
open Lean Meta Elab Command

#eval show CommandElabM Unit from do
  let env ← getEnv
  let allowed : List Name := [``propext, ``Classical.choice, ``Quot.sound]
  let mut bad : Array (Name × Array Name) := #[]
  let mut count : Nat := 0
  for (n, ci) in env.constants.toList do
    if let some idx := env.getModuleIdxFor? n then
      let modName := env.header.moduleNames[idx.toNat]!
      if (`DFL).isPrefixOf modName then
        if ci.isTheorem then
          count := count + 1
          let axs ← liftCoreM <| Lean.collectAxioms n
          if axs.any (fun a => !(allowed.contains a)) then
            bad := bad.push (n, axs)
  logInfo m!"checked {count} theorems in DFL modules; with non-standard axioms: {bad.size}"
  for (n, axs) in bad do
    logInfo m!"{n}: {axs}"
