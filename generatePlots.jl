#author: Norbert Jaśniewicz

import DataFrames, CSV
using Plots


function main()
  plots_names = ["LU", "Gauss", "LU pivot", "Gauss pivot"]

  t = CSV.read("experimentsTime.csv", DataFrames.DataFrame)
  p = Plots.plot(legend=:topleft, xlabel="Size", ylabel="Time[s]", title="Time complexity", dpi=600)
  for i in (2:size(t)[2])
    Plots.plot!(t[:,1], t[:,i], label=plots_names[i - 1])
  end  
  savefig(p, "time_complexity")

  t = CSV.read("experimentsMemory.csv", DataFrames.DataFrame)
  p = Plots.plot(legend=:topleft, xlabel="Size", ylabel="Allocated memory[B]", title="Space complexity", dpi=600)
  for i in (2:size(t)[2])
    Plots.plot!(t[:,1], t[:,i], label=plots_names[i - 1])
  end  
  savefig(p, "space_complexity")

  t = CSV.read("experimentsLoad.csv", DataFrames.DataFrame)
  p = Plots.plot(legend=:topleft, xlabel="Size", ylabel="Allocated memory[B]", title="Space complexity of loaded matrix", dpi=600)
  for i in (2:size(t)[2])
    Plots.plot!(t[:,1], t[:,i], label="SparseMatrix")
  end  
  savefig(p, "load_complexity")
end

main()