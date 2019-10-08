
evaluation.metric = function( Net.true, Net.pred){

  t = 1e-5

  Net.pred= (abs(Net.pred) > t)*upper.tri(Net.pred)
  Net.true =(abs(Net.true) > t)*upper.tri(Net.true)

  TP = sum((Net.true!=0)*(Net.pred!=0))
  FP = sum((Net.true==0)*(Net.pred!=0))
  TN = sum((Net.true==0)*(Net.pred==0))
  FN = sum((Net.true!=0)*(Net.pred==0))

  result = list( TP = TP,  FP = FP , TN = TN, FN = FN)
}
