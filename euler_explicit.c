static char help[] = "Euler explicit \n\n";

#include <petscmat.h>
#include <stdio.h>

int main(int argc,char **args)
{  
  Mat            A; 
  Vec            u,uold,f;
  PetscReal      error=1000,tol=1e-8;/* norm of solution error */
  PetscReal      rho,c,k;
  PetscReal      l,t,x,dt,dx;
  PetscReal      p,a;
  PetscReal      val,zero = 0.0;
  PetscErrorCode ierr;
  PetscInt       i,n=10,nn,col[3];
  PetscInt       iter,rstart,rend,nlocal;
  PetscScalar    value[3],normU,normUold;
  
  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"\nmesh size n=%d \n",n);
  
  rho = 1.0;
  c = 1.0;
  k = 1.0;
  l = 1.0;
  x = 0.0;
  t = 0.0;
  dx = l/n;
  dt = 0.005;
  iter = 0;
  
  nn = n+1;
  ierr = VecCreate(PETSC_COMM_WORLD,&u);CHKERRQ(ierr);
  ierr = VecSetSizes(u,PETSC_DECIDE,nn);CHKERRQ(ierr);
  ierr = VecSetFromOptions(u);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&uold);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&f);CHKERRQ(ierr);
  
  ierr = VecGetOwnershipRange(u,&rstart,&rend);CHKERRQ(ierr);
  ierr = VecGetLocalSize(u,&nlocal);CHKERRQ(ierr);
  
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,nlocal,nlocal,nn,nn);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);
  
  /* initialize vector uold,f,u */
  for(i=rstart;i<rend;i++){
	  x = i*dx;
	  val = PetscSinReal(PETSC_PI*x)*dt/(rho*c);
	  VecSetValues(f,1,&i,&val,INSERT_VALUES);
	  val = PetscExpReal(x);
	  VecSetValues(uold,1,&i,&val,INSERT_VALUES);
  }
  ierr = VecSet(u,zero);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(u);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(u);CHKERRQ(ierr);
//ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(uold);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(uold);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(f);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(f);CHKERRQ(ierr);
  
  //set boundary conditions
   VecSetValue(uold,0,0.0,INSERT_VALUES);// BC
   VecSetValue(uold,n,0.0,INSERT_VALUES);// BC
   ierr = VecAssemblyBegin(uold);CHKERRQ(ierr);
   ierr = VecAssemblyEnd(uold);CHKERRQ(ierr);
ierr = VecView(uold,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); 
   VecSetValue(f,0,0.0,INSERT_VALUES);// BC
   VecSetValue(f,n,0.0,INSERT_VALUES);// BC
   ierr = VecAssemblyBegin(f);CHKERRQ(ierr);
   ierr = VecAssemblyEnd(f);CHKERRQ(ierr); 
ierr = VecView(f,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); 
  
  /*initialize matrix A */
  p = k/(rho*c);
  a = (p*dt)/(dx*dx);
  if (!rstart) 
  {
    rstart = 1;
    i      = 0; col[0] = 0; col[1] = 1; value[0] = 1; value[1] = 0;
    ierr   = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }
  if (rend == nn) 
  {
    rend = nn-1;
    i    = nn-1; col[0] = nn-2; col[1] = nn-1; value[0] = 0; value[1] = 1;
    ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }
  /* Set entries corresponding to the mesh interior */
  value[0] = a; value[1] = 1-2*a; value[2] = a;
  for (i=rstart; i<rend; i++) 
  {
    col[0] = i-1; col[1] = i; col[2] = i+1;
    ierr   = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  
      
  while(error>tol){
	  error = 0.0;
	  ierr = MatMult(A,uold,u);CHKERRQ(ierr);
	  ierr = VecAYPX(u,1,f);CHKERRQ(ierr);
//ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	  //VecSetValue(u,0,0.0,INSERT_VALUES);// BC
	  //VecSetValue(u,n,0.0,INSERT_VALUES);// BC
	  //ierr = VecAssemblyBegin(u);CHKERRQ(ierr);
          //ierr = VecAssemblyEnd(u);CHKERRQ(ierr);
	  ierr = VecNorm(uold,NORM_2,&normUold);CHKERRQ(ierr);
	  ierr = VecNorm(u,NORM_2,&normU);CHKERRQ(ierr);
	  error=PetscAbsReal((normU-normUold)/normUold);
	  ierr = VecAYPX(uold,0,u);CHKERRQ(ierr);
	  t = t + dt;
	  iter = iter + 1;
  }
  
  
  ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"iter=%d, error=%.6g \n",iter,error);
  
  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&uold);CHKERRQ(ierr);
  ierr = VecDestroy(&f);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  
  
  ierr = PetscFinalize();
  return ierr;
  
}
