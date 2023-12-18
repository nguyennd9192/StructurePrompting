
#!/bin/sh
code_dir="code_dir"
i=1
while [ ! -f ./USPEX_IS_DONE ] ; do
   date >> log
   USPEX -r >> log
   python $code_dir/collect_structures.py example/Sm1Fe12
   python $code_dir/descriptor.py example/Sm1Fe12
   python $code_dir/process_features.py example/Sm1Fe12
   python $code_dir/collect_enthalpy_force.py example/Sm1Fe12
   python $code_dir/seeding.py example/Sm1Fe12
   sleep 500
done
