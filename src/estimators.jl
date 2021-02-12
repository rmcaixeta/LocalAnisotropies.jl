
@estimsolver LocalKriging begin
  @param variogram = (:X => ExponentialVariogram())
  @param mean = nothing #0.0 or [1.0,0.8,....]
  @param method = :MovingWindows
  @param localpars = nothing
  @param localparshd = nothing
  @param minneighbors = 1
  @param maxneighbors = 20
  @param neighborhood = nothing
  @param distance = Euclidean()
end
