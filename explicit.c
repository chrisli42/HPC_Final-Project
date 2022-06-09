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
  
/* HDF5: Initialization */
	hid_t        file_id, dataset_id1, dataset_id2, group_id, dataspace_id1, dataspace_id2;  /* identifiers */
	hsize_t      dim1[1], dim2[1];
	herr_t       status;
	double* vec1 = (double*)malloc(N * sizeof(double));
	double* vec2 = (double*)malloc(3 * sizeof(double));
	free(vec1);
	free(vec2);

	dim1[0] = N;
	dim2[0] = 3;





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

  

bool restart = false;
	if (restart == false)
	{
		/* HDF5: Create a new file to store datasets. */
		file_id = H5Fcreate(FILE2, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

		/* HDF5: data1-output velocity
				 data2-parameters:t, dt, N */
		dataspace_id1 = H5Screate_simple(1, dim1, NULL);
		dataspace_id2 = H5Screate_simple(1, dim2, NULL);

		/* Create two groups to store initial values and output values in the file. */
		group_id = H5Gcreate2(file_id, "/output", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		/* Create the datasets. */
		dataset_id1 = H5Dcreate2(file_id, "/output/uout", H5T_IEEE_F64BE, dataspace_id1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		dataset_id2 = H5Dcreate2(file_id, "/output/parameters", H5T_IEEE_F64BE, dataspace_id2, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		para[0] = t;
		para[1] = dt;
		para[2] = n;

		// Set values to u
		ierr = VecGetOwnershipRange(u, &rstart, &rend); CHKERRQ(ierr);
		for (int i = rstart; i < rend; i++) {
			temp = exp(((double)i - 0.5) * dx);
			ierr = VecSetValue(u, i, temp, INSERT_VALUES); CHKERRQ(ierr);
		}
		ierr = VecAssemblyBegin(u); CHKERRQ(ierr);
		ierr = VecAssemblyEnd(u); CHKERRQ(ierr);
		//ierr = VecView(u, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	}

	else
	{
		/* HDF5: Open an existing file which stores datasets. */
		file_id = H5Fopen(FILE2, H5F_ACC_RDWR, H5P_DEFAULT);

		/* HDF5: data1-output velocity
				 data2-parameters:t, dt, N */

				 /* Create two groups to store initial values and output values in the file. */
		group_id = H5Gopen(file_id, "/output", H5P_DEFAULT);

		/* Create the datasets. */
		dataset_id1 = H5Dopen(file_id, "/output/uout", H5P_DEFAULT);
		dataset_id2 = H5Dopen(file_id, "/output/parameters", H5P_DEFAULT);

		/* Write data from hdf5 to arrays */
		status = H5Dread(dataset_id1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, uout);
		status = H5Dread(dataset_id2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, para);

		t = para[0];
		dt = para[1];
		n = para[2];

		// Set values to u-from uout in hdf5
		ierr = VecGetOwnershipRange(u, &rstart, &rend); CHKERRQ(ierr);
		for (int i = rstart; i < rend; i++) {
			temp = uout[i];
			ierr = VecSetValue(u, i, temp, INSERT_VALUES); CHKERRQ(ierr);
		}
		ierr = VecAssemblyBegin(u); CHKERRQ(ierr);
		ierr = VecAssemblyEnd(u); CHKERRQ(ierr);
		//ierr = VecView(u, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	}





      
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

  if ((iter + 1) % 10 == 0)  // HDf5 record every 10 step
		{
			for (int i = 0; i < N; i++)
			{
				num[i] = i;
			}
			ierr = VecGetValues(u, N, num, uout); CHKERRQ(ierr);
			status = H5Dwrite(dataset_id1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, uout);
			status = H5Dwrite(dataset_id2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, para);
		}




  }
  
  
  ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"iter=%d, error=%.6g \n",iter,error);
  
  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&uold);CHKERRQ(ierr);
  ierr = VecDestroy(&f);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  
  
  ierr = PetscFinalize();
  return ierr;
  


  status = H5Sclose(dataspace_id1);
	status = H5Sclose(dataspace_id2);

	/* Close the first dataset. */
	status = H5Dclose(dataset_id1);
	status = H5Dclose(dataset_id2);
	/* Close the group. */
	status = H5Gclose(group_id);
	status = H5Fclose(file_id);


}
