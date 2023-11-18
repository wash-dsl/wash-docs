rm doc_out/sphinx/*.html 
rm api/*
doxygen Doxyfile.in && sphinx-build -b html -Dbreathe_projects.my_project=doc_out/xml . doc_out/sphinx/ && cp -a doc_out/sphinx/. ../docs
