# simple tests of TaskInterface

import TaskInterface
# test NOBAT
xxx=TaskInterface.TaskInterface("NOBAT","NOBAT.TDF")
xxx.Auserid=100
xxx.setValue("DETIME",1.0)
xxx.Input()
xxx.Check()
xxx.Go()

# test IMEAN
xxx=TaskInterface.TaskInterface("IMEAN","IMEAN.TDF")
xxx.Auserid=102
xxx.setValue("DOHIST",1.0)
xxx.setValue("USERID",102.0)
xxx.setValue("INNAME","L1_AB_TOT")
xxx.setValue("INCLASS","SUBIM")
xxx.setValue("INDISK",1.0)
xxx.setValue("INSEQ",10.0)
xxx.setValue("OUTFILE","python.list")
xxx.Input()
xxx.Check()
xxx.Go()

