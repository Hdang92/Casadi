nx =nlp.x;
nf=nlp.f;
ng = nlp.g ;

jac_f_fun = Function('jac_f_fun',{nlp.x, nlp.p},{jacobian(nlp.f,nlp.x)}); 
jac_g_fun = Function('jac_g_fun',{nlp.x, nlp.p},{jacobian(nlp.g,nlp.x)}); 
jac_f_opt = jac_f_fun(sol.x, m);
jac_g_opt = jac_g_fun(sol.x, m); 
M = [jac_f_opt'*jac_f_opt jac_g_opt'; jac_g_opt zeros(size(nlp.g,1))];

Minv = solve(M, [eye(size(nlp.x,1));zeros(size(jac_g_opt))]); 
Bock2007b 
Chapter 4 