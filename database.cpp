#include "database.h"

// initialize it for singleton. otherwise will get LNK2001 error
database* database::instance = 0;
int database::refCount = 0;

database::database(void)
{
	textureimage=NULL;
	themesh=NULL;
	cubic_bezier=NULL;
}


database::~database(void)
{
	//if(&seeds_golobal_coord)
	//	seeds_golobal_coord.clear();
	clear();
}

database* database::getInstance()
{
	if(refCount == 0){
		instance = new database();
	    refCount++; //increase reference counter
	}
	return instance;
}

void database::Release()
{
	if( --refCount < 1 )
	{
		Destroy();
	}
}

void database::Destroy()
{
	if( instance != 0 )
	{
		delete( instance );
		instance = 0;
	}
}

void database::PN_TRIANGLE(int level_num){
	cubic_bezier =new TriMesh;
	cubic_bezier->request_vertex_normals();
	TriMesh::Point the;
	TriMesh::Normal the_norm;
	TriMesh::VertexHandle vh_the;
	std::vector<TriMesh::VertexHandle> vh;
	std::vector<TriMesh::VertexHandle> temp_fh;
	
	TriMesh::FaceIter f_it=themesh->faces_begin(), f_end=themesh->faces_end();
	for ( ; f_it!=f_end; ++f_it)	{

		vh.clear();

		TriMesh::FaceVertexIter fv_it=themesh->fv_iter(f_it.handle());
		TriMesh::Point b300=themesh->point(fv_it.handle());
		TriMesh::Normal N1=themesh->normal(fv_it.handle());
		vh_the=cubic_bezier->add_vertex(b300);
		cubic_bezier->set_normal(vh_the,N1);
		vh.push_back(vh_the);
		++fv_it;
		TriMesh::Point b030=themesh->point(fv_it.handle());
		TriMesh::Normal N2=themesh->normal(fv_it.handle());
		vh_the=cubic_bezier->add_vertex(b030);
		cubic_bezier->set_normal(vh_the,N2);
		vh.push_back(vh_the);
		++fv_it;
		TriMesh::Point b003=themesh->point(fv_it.handle());
		TriMesh::Normal N3=themesh->normal(fv_it.handle());
		vh_the=cubic_bezier->add_vertex(b003);
		cubic_bezier->set_normal(vh_the,N3);
		vh.push_back(vh_the);

		TriMesh::Point P1=b300; TriMesh::Point P2=b030; TriMesh::Point P3=b003;

		double w12=DOTPROD3(P2-P1, N1);
		double w21=DOTPROD3(P1-P2, N2);
		double w13=DOTPROD3(P3-P1, N1);
		double w31=DOTPROD3(P1-P3, N3);
		double w23=DOTPROD3(P3-P2, N2);
		double w32=DOTPROD3(P2-P3, N3);

		TriMesh::Point b210=(2.0*P1+P2-w12*N1)/3.0;
		TriMesh::Point b120=(2.0*P2+P1-w21*N2)/3.0;
		TriMesh::Point b021=(2.0*P2+P3-w23*N2)/3.0;
		TriMesh::Point b012=(2.0*P3+P2-w32*N3)/3.0;
		TriMesh::Point b102=(2.0*P3+P1-w31*N3)/3.0;
		TriMesh::Point b201=(2.0*P1+P3-w13*N1)/3.0;
		TriMesh::Point E=(b210+b120+b012+b021+b102+b201)/6.0;
		TriMesh::Point V=(P1+P2+P3)/3.0;
		TriMesh::Point b111=E+(E-V)/2.0;
		
		double v12=2.0*DOTPROD3(P2-P1, N1+N2)/DOTPROD3(P2-P1, P2-P1);
		double v23=2.0*DOTPROD3(P3-P2, N2+N3)/DOTPROD3(P3-P2, P3-P2);
		double v31=2.0*DOTPROD3(P1-P3, N3+N1)/DOTPROD3(P1-P3, P1-P3);

		TriMesh::Normal n110=N1+N2-v12*(P2-P1); GeometryTools::normalize(n110);
		TriMesh::Normal n011=N2+N3-v23*(P3-P2); GeometryTools::normalize(n011);
		TriMesh::Normal n101=N3+N1-v31*(P1-P3); GeometryTools::normalize(n101);

		TriMesh::Normal n210=N1*2.0*2.0/3.0/3.0+N2/3.0/3.0+n110*2.0/3.0/3.0;
		GeometryTools::normalize(n210);
		TriMesh::Normal n120=N2*2.0*2.0/3.0/3.0+N1/3.0/3.0+n110*2.0/3.0/3.0;
		GeometryTools::normalize(n120);

		TriMesh::Normal n021=N2*2.0*2.0/3.0/3.0+N3/3.0/3.0+n011*2.0/3.0/3.0;
		GeometryTools::normalize(n021);
		TriMesh::Normal n012=N3*2.0*2.0/3.0/3.0+N2/3.0/3.0+n011*2.0/3.0/3.0;
		GeometryTools::normalize(n012);

		TriMesh::Normal n201=N1*2.0*2.0/3.0/3.0+N3/3.0/3.0+n101*2.0/3.0/3.0;
		GeometryTools::normalize(n201);
		TriMesh::Normal n102=N3*2.0*2.0/3.0/3.0+N1/3.0/3.0+n101*2.0/3.0/3.0;
		GeometryTools::normalize(n102);

		TriMesh::Normal n111=N1/3.0/3.0+N2/3.0/3.0+N3/3.0/3.0+n101/3.0/3.0+n110/3.0/3.0+n011/3.0/3.0;
		GeometryTools::normalize(n111);

		
		vh_the=cubic_bezier->add_vertex(b120);
		cubic_bezier->set_normal(vh_the,n120);
		vh.push_back(vh_the);
		vh_the=cubic_bezier->add_vertex(b021);
		cubic_bezier->set_normal(vh_the,n021);
		vh.push_back(vh_the);

		vh_the=cubic_bezier->add_vertex(b210);
		cubic_bezier->set_normal(vh_the,n210);
		vh.push_back(vh_the);
		vh_the=cubic_bezier->add_vertex(b111);
		cubic_bezier->set_normal(vh_the,n111);
		vh.push_back(vh_the);
		vh_the=cubic_bezier->add_vertex(b012);
		cubic_bezier->set_normal(vh_the,n012);
		vh.push_back(vh_the);

		vh_the=cubic_bezier->add_vertex(b201);
		cubic_bezier->set_normal(vh_the,n201);
		vh.push_back(vh_the);
		vh_the=cubic_bezier->add_vertex(b102);
		cubic_bezier->set_normal(vh_the,n102);
		vh.push_back(vh_the);


		//row 1
		temp_fh.clear();
		temp_fh.push_back(vh[3]); temp_fh.push_back(vh[1]); temp_fh.push_back(vh[4]);
		cubic_bezier->add_face(temp_fh);
		temp_fh.clear();
		//row 2
		temp_fh.push_back(vh[5]); temp_fh.push_back(vh[3]); temp_fh.push_back(vh[6]);
		cubic_bezier->add_face(temp_fh);
		temp_fh.clear();
		temp_fh.push_back(vh[6]); temp_fh.push_back(vh[3]); temp_fh.push_back(vh[4]);
		cubic_bezier->add_face(temp_fh);
		temp_fh.clear();
		temp_fh.push_back(vh[6]); temp_fh.push_back(vh[4]); temp_fh.push_back(vh[7]);
		cubic_bezier->add_face(temp_fh);
		temp_fh.clear();
		//row 3
		temp_fh.push_back(vh[0]); temp_fh.push_back(vh[5]); temp_fh.push_back(vh[8]);
		cubic_bezier->add_face(temp_fh);
		temp_fh.clear();
		temp_fh.push_back(vh[8]); temp_fh.push_back(vh[5]); temp_fh.push_back(vh[6]);
		cubic_bezier->add_face(temp_fh);
		temp_fh.clear();
		temp_fh.push_back(vh[8]); temp_fh.push_back(vh[6]); temp_fh.push_back(vh[9]);
		cubic_bezier->add_face(temp_fh);
		temp_fh.clear();
		temp_fh.push_back(vh[9]); temp_fh.push_back(vh[6]); temp_fh.push_back(vh[7]);
		cubic_bezier->add_face(temp_fh);
		temp_fh.clear();
		temp_fh.push_back(vh[9]); temp_fh.push_back(vh[7]); temp_fh.push_back(vh[2]);
		cubic_bezier->add_face(temp_fh);
		temp_fh.clear();
	}
	vh.clear();


	for(int i=2;i<=level_num;i++){
		TriMesh *temp_mesh=new TriMesh;
		temp_mesh->request_vertex_normals();

		TriMesh::FaceIter f_it=cubic_bezier->faces_begin(), f_end=cubic_bezier->faces_end();
		for ( ; f_it!=f_end; ++f_it)	{

			vh.clear();

			TriMesh::FaceVertexIter fv_it=cubic_bezier->fv_iter(f_it.handle());
			TriMesh::Point b300=cubic_bezier->point(fv_it.handle());
			TriMesh::Normal N1=cubic_bezier->normal(fv_it.handle());
			vh_the=temp_mesh->add_vertex(b300);
			temp_mesh->set_normal(vh_the,N1);
			vh.push_back(vh_the);
			++fv_it;
			TriMesh::Point b030=cubic_bezier->point(fv_it.handle());
			TriMesh::Normal N2=cubic_bezier->normal(fv_it.handle());
			vh_the=temp_mesh->add_vertex(b030);
			temp_mesh->set_normal(vh_the,N2);
			vh.push_back(vh_the);
			++fv_it;
			TriMesh::Point b003=cubic_bezier->point(fv_it.handle());
			TriMesh::Normal N3=cubic_bezier->normal(fv_it.handle());
			vh_the=temp_mesh->add_vertex(b003);
			temp_mesh->set_normal(vh_the,N3);
			vh.push_back(vh_the);

			temp_mesh->add_face(vh);
		}

		cubic_bezier->clean();
        
		f_it=temp_mesh->faces_begin(); f_end=temp_mesh->faces_end();
		for ( ; f_it!=f_end; ++f_it)	{

			vh.clear();

			TriMesh::FaceVertexIter fv_it=temp_mesh->fv_iter(f_it.handle());
			TriMesh::Point b300=temp_mesh->point(fv_it.handle());
			TriMesh::Normal N1=temp_mesh->normal(fv_it.handle());
			vh_the=cubic_bezier->add_vertex(b300);
			cubic_bezier->set_normal(vh_the,N1);
			vh.push_back(vh_the);
			++fv_it;
			TriMesh::Point b030=temp_mesh->point(fv_it.handle());
			TriMesh::Normal N2=temp_mesh->normal(fv_it.handle());
			vh_the=cubic_bezier->add_vertex(b030);
			cubic_bezier->set_normal(vh_the,N2);
			vh.push_back(vh_the);
			++fv_it;
			TriMesh::Point b003=temp_mesh->point(fv_it.handle());
			TriMesh::Normal N3=temp_mesh->normal(fv_it.handle());
			vh_the=cubic_bezier->add_vertex(b003);
			cubic_bezier->set_normal(vh_the,N3);
			vh.push_back(vh_the);

			TriMesh::Point P1=b300; TriMesh::Point P2=b030; TriMesh::Point P3=b003;

			double w12=DOTPROD3(P2-P1, N1);
			double w21=DOTPROD3(P1-P2, N2);
			double w13=DOTPROD3(P3-P1, N1);
			double w31=DOTPROD3(P1-P3, N3);
			double w23=DOTPROD3(P3-P2, N2);
			double w32=DOTPROD3(P2-P3, N3);

			TriMesh::Point b210=(2.0*P1+P2-w12*N1)/3.0;
			TriMesh::Point b120=(2.0*P2+P1-w21*N2)/3.0;
			TriMesh::Point b021=(2.0*P2+P3-w23*N2)/3.0;
			TriMesh::Point b012=(2.0*P3+P2-w32*N3)/3.0;
			TriMesh::Point b102=(2.0*P3+P1-w31*N3)/3.0;
			TriMesh::Point b201=(2.0*P1+P3-w13*N1)/3.0;
			TriMesh::Point E=(b210+b120+b012+b021+b102+b201)/6.0;
			TriMesh::Point V=(P1+P2+P3)/3.0;
			TriMesh::Point b111=E+(E-V)/2.0;

			double v12=2.0*DOTPROD3(P2-P1, N1+N2)/DOTPROD3(P2-P1, P2-P1);
			double v23=2.0*DOTPROD3(P3-P2, N2+N3)/DOTPROD3(P3-P2, P3-P2);
			double v31=2.0*DOTPROD3(P1-P3, N3+N1)/DOTPROD3(P1-P3, P1-P3);

			TriMesh::Normal n110=N1+N2-v12*(P2-P1); GeometryTools::normalize(n110);
			TriMesh::Normal n011=N2+N3-v23*(P3-P2); GeometryTools::normalize(n011);
			TriMesh::Normal n101=N3+N1-v31*(P1-P3); GeometryTools::normalize(n101);

			TriMesh::Normal n210=N1*2.0*2.0/3.0/3.0+N2/3.0/3.0+n110*2.0/3.0/3.0;
			GeometryTools::normalize(n210);
			TriMesh::Normal n120=N2*2.0*2.0/3.0/3.0+N1/3.0/3.0+n110*2.0/3.0/3.0;
			GeometryTools::normalize(n120);

			TriMesh::Normal n021=N2*2.0*2.0/3.0/3.0+N3/3.0/3.0+n011*2.0/3.0/3.0;
			GeometryTools::normalize(n021);
			TriMesh::Normal n012=N3*2.0*2.0/3.0/3.0+N2/3.0/3.0+n011*2.0/3.0/3.0;
			GeometryTools::normalize(n012);

			TriMesh::Normal n201=N1*2.0*2.0/3.0/3.0+N3/3.0/3.0+n101*2.0/3.0/3.0;
			GeometryTools::normalize(n201);
			TriMesh::Normal n102=N3*2.0*2.0/3.0/3.0+N1/3.0/3.0+n101*2.0/3.0/3.0;
			GeometryTools::normalize(n102);

			TriMesh::Normal n111=N1/3.0/3.0+N2/3.0/3.0+N3/3.0/3.0+n101/3.0/3.0+n110/3.0/3.0+n011/3.0/3.0;
			GeometryTools::normalize(n111);


			vh_the=cubic_bezier->add_vertex(b120);
			cubic_bezier->set_normal(vh_the,n120);
			vh.push_back(vh_the);
			vh_the=cubic_bezier->add_vertex(b021);
			cubic_bezier->set_normal(vh_the,n021);
			vh.push_back(vh_the);

			vh_the=cubic_bezier->add_vertex(b210);
			cubic_bezier->set_normal(vh_the,n210);
			vh.push_back(vh_the);
			vh_the=cubic_bezier->add_vertex(b111);
			cubic_bezier->set_normal(vh_the,n111);
			vh.push_back(vh_the);
			vh_the=cubic_bezier->add_vertex(b012);
			cubic_bezier->set_normal(vh_the,n012);
			vh.push_back(vh_the);

			vh_the=cubic_bezier->add_vertex(b201);
			cubic_bezier->set_normal(vh_the,n201);
			vh.push_back(vh_the);
			vh_the=cubic_bezier->add_vertex(b102);
			cubic_bezier->set_normal(vh_the,n102);
			vh.push_back(vh_the);


			//row 1
			temp_fh.clear();
			temp_fh.push_back(vh[3]); temp_fh.push_back(vh[1]); temp_fh.push_back(vh[4]);
			cubic_bezier->add_face(temp_fh);
			temp_fh.clear();
			//row 2
			temp_fh.push_back(vh[5]); temp_fh.push_back(vh[3]); temp_fh.push_back(vh[6]);
			cubic_bezier->add_face(temp_fh);
			temp_fh.clear();
			temp_fh.push_back(vh[6]); temp_fh.push_back(vh[3]); temp_fh.push_back(vh[4]);
			cubic_bezier->add_face(temp_fh);
			temp_fh.clear();
			temp_fh.push_back(vh[6]); temp_fh.push_back(vh[4]); temp_fh.push_back(vh[7]);
			cubic_bezier->add_face(temp_fh);
			temp_fh.clear();
			//row 3
			temp_fh.push_back(vh[0]); temp_fh.push_back(vh[5]); temp_fh.push_back(vh[8]);
			cubic_bezier->add_face(temp_fh);
			temp_fh.clear();
			temp_fh.push_back(vh[8]); temp_fh.push_back(vh[5]); temp_fh.push_back(vh[6]);
			cubic_bezier->add_face(temp_fh);
			temp_fh.clear();
			temp_fh.push_back(vh[8]); temp_fh.push_back(vh[6]); temp_fh.push_back(vh[9]);
			cubic_bezier->add_face(temp_fh);
			temp_fh.clear();
			temp_fh.push_back(vh[9]); temp_fh.push_back(vh[6]); temp_fh.push_back(vh[7]);
			cubic_bezier->add_face(temp_fh);
			temp_fh.clear();
			temp_fh.push_back(vh[9]); temp_fh.push_back(vh[7]); temp_fh.push_back(vh[2]);
			cubic_bezier->add_face(temp_fh);
			temp_fh.clear();
		}
		vh.clear();

		temp_mesh->clear();
		delete []temp_mesh;
	}
}

void database::exact_pn_triangle(int level_num){
	cubic_bezier =new TriMesh;
	cubic_bezier->request_vertex_normals();
	TriMesh::Point the;
	TriMesh::Normal the_norm;
	TriMesh::VertexHandle vh_the;
	std::vector<TriMesh::VertexHandle> vh;
	std::vector<TriMesh::VertexHandle> temp_fh;

	TriMesh::FaceIter f_it=themesh->faces_begin(), f_end=themesh->faces_end();
	for ( ; f_it!=f_end; ++f_it)	{
		TriMesh::FaceVertexIter fv_it=themesh->fv_iter(f_it.handle());
		TriMesh::Point b300=themesh->point(fv_it.handle());
		TriMesh::Normal N1=themesh->normal(fv_it.handle());
		++fv_it;
		TriMesh::Point b030=themesh->point(fv_it.handle());
		TriMesh::Normal N2=themesh->normal(fv_it.handle());
		++fv_it;
		TriMesh::Point b003=themesh->point(fv_it.handle());
		TriMesh::Normal N3=themesh->normal(fv_it.handle());

		TriMesh::Point P1=b300; TriMesh::Point P2=b030; TriMesh::Point P3=b003;

		double w12=DOTPROD3(P2-P1, N1);
		double w21=DOTPROD3(P1-P2, N2);
		double w13=DOTPROD3(P3-P1, N1);
		double w31=DOTPROD3(P1-P3, N3);
		double w23=DOTPROD3(P3-P2, N2);
		double w32=DOTPROD3(P2-P3, N3);

		TriMesh::Point b210=(2.0*P1+P2-w12*N1)/3.0;
		TriMesh::Point b120=(2.0*P2+P1-w21*N2)/3.0;
		TriMesh::Point b021=(2.0*P2+P3-w23*N2)/3.0;
		TriMesh::Point b012=(2.0*P3+P2-w32*N3)/3.0;
		TriMesh::Point b102=(2.0*P3+P1-w31*N3)/3.0;
		TriMesh::Point b201=(2.0*P1+P3-w13*N1)/3.0;
		TriMesh::Point E=(b210+b120+b012+b021+b102+b201)/6.0;
		TriMesh::Point V=(P1+P2+P3)/3.0;
		TriMesh::Point b111=E+(E-V)/2.0;

		double v12=2.0*DOTPROD3(P2-P1, N1+N2)/DOTPROD3(P2-P1, P2-P1);
		double v23=2.0*DOTPROD3(P3-P2, N2+N3)/DOTPROD3(P3-P2, P3-P2);
		double v31=2.0*DOTPROD3(P1-P3, N3+N1)/DOTPROD3(P1-P3, P1-P3);

		TriMesh::Normal n110=N1+N2-v12*(P2-P1); GeometryTools::normalize(n110);
		TriMesh::Normal n011=N2+N3-v23*(P3-P2); GeometryTools::normalize(n011);
		TriMesh::Normal n101=N3+N1-v31*(P1-P3); GeometryTools::normalize(n101);

		double inc=1.0/(level_num+1);
		double w=0.0; double u=0.0; double v=0.0;
		//generate each vertex
		for (int i=level_num+1; i>=0;i--){
			w=(i+0.0)/(level_num+1+0.0);
			v=0.0;
			u=1.0-w-v;
			while(u>-0.00000000001){
				TriMesh::Point temp_point=b300*w*w*w+b030*u*u*u+b003*v*v*v
					+3.0*b120*w*u*u+3.0*b210*w*w*u
					+3.0*b201*w*w*v+3.0*b102*w*v*v
					+3.0*b012*u*v*v+3.0*b021*u*u*v
					+6.0*b111*w*u*v;
				TriMesh::Normal temp_normal=N1*w*w+N2*u*u+N3*v*v+n110*w*u+n011*u*v+n101*w*v;
				GeometryTools::normalize(temp_normal);
				vh_the=cubic_bezier->add_vertex(temp_point);
				cubic_bezier->set_normal(vh_the,temp_normal);
				vh.push_back(vh_the);

				v+=inc;
				u-=inc;
			}
		}
		//generate each triangle line by line
		for (int i=0; i<=level_num;i++){
			int up_index=i*(i+1)/2;
			int down_index=(i+1)*(i+2)/2;

			for(int j=0; j<i+1;j++){
				temp_fh.push_back(vh[up_index]); temp_fh.push_back(vh[down_index]); temp_fh.push_back(vh[down_index+1]);
				cubic_bezier->add_face(temp_fh);
				temp_fh.clear();
				up_index++;
				down_index++;
			}

			up_index=i*(i+1)/2;
			down_index=(i+1)*(i+2)/2;
			for(int j=0; j<i;j++){
				temp_fh.push_back(vh[up_index]); temp_fh.push_back(vh[down_index+1]); temp_fh.push_back(vh[up_index+1]); 
				cubic_bezier->add_face(temp_fh);
				temp_fh.clear();
				up_index++;
				down_index++;
			}
		}//generate each triangle line by line
		vh.clear();
	}// for each triangle of the original triangle in trimesh
}

void database::exact_torus(int level_num){
	cubic_bezier =new TriMesh;
	cubic_bezier->request_vertex_normals();
	TriMesh::Point the;
	TriMesh::Normal the_norm;
	TriMesh::VertexHandle vh_the;
	std::vector<TriMesh::VertexHandle> vh;
	std::vector<TriMesh::VertexHandle> temp_fh;

	double u1=(int)(2/12)*PI/12.5; double v1=PI/3.0;
	double u2=(int)(10/12)*PI/12.5; double v2=5.0*PI/3.0;
	double u3=(int)(24/12)*PI/12.5; double v3=0.0;

	double u=0.0;
	double v=0.0;
	for (int i=0; i<level_num+1;i++){
		u=u3-i*(u3-u1)/level_num;
		v=v3+(v1-v3)*i/level_num;
		for(int j=0; j<level_num+1;j++){
			TriMesh::Point temp_point;
			temp_point[0]=(3.0+cos(v))*cos(u);
			temp_point[1]=(3.0+cos(v))*sin(u);
			temp_point[2]=sin(v);

			TriMesh::Normal d_u,d_v;
			d_u[0]=-(3.0+cos(v))*sin(u);
			d_u[1]=(3.0+cos(v))*cos(u);
			d_u[2]=0.0;
			d_v[0]=-cos(u)*sin(v);
			d_v[1]=-sin(u)*sin(v);
			d_v[2]=cos(v);
			TriMesh::Normal temp_normal;
			CROSSPROD3(temp_normal, d_u,d_v);
			GeometryTools::normalize(temp_normal);
			vh_the=cubic_bezier->add_vertex(temp_point);
			cubic_bezier->set_normal(vh_the,temp_normal);
			vh.push_back(vh_the);

			v-=(v1+2.0*PI-v2)/level_num;
			if(v<0.0)
				v+=2.0*PI;
		}
	}
	//generate each triangle line by line
	for (int i=0; i<level_num;i++){
		int up_index=i*(level_num+1);
		int down_index=(i+1)*(level_num+1);

		for(int j=0; j<level_num;j++){
			temp_fh.push_back(vh[up_index]); temp_fh.push_back(vh[down_index]); temp_fh.push_back(vh[down_index+1]);
			cubic_bezier->add_face(temp_fh);
			temp_fh.clear();
			temp_fh.push_back(vh[up_index+1]);  temp_fh.push_back(vh[up_index]);   temp_fh.push_back(vh[down_index+1]); 
			cubic_bezier->add_face(temp_fh);
			temp_fh.clear();
			up_index++;
			down_index++;
		}
	}//generate each triangle line by line
	vh.clear();

		//double u1=(int)(2/12)*PI/12.5; double v1=PI/3.0;
		//double u2=(int)(10/12)*PI/12.5; double v2=5.0*PI/3.0;
		//double u3=(int)(24/12)*PI/12.5; double v3=0.0;
		////startu=(int)((startindex)/12)*PI/12.5;
		////startv=startindex%12;
		////startv=startv*PI/6.0;
		////generate each vertex
		//double u=0.0;
		//double v=0.0;
		//for (int i=0; i<level_num+1;i++){
		//	//u=u3+i*(u3-u1)/level_num;
		//	//v=v3+(v1+2.0*PI-v2)*i/level_num;
		//	u=u3-i*(u3-u1)/level_num;
		//	v=v3+(v1-v3)*i/level_num;
		//	for(int j=0; j<i+1;j++){
		//		TriMesh::Point temp_point;
		//		temp_point[0]=(3.0+cos(v))*cos(u);
		//		temp_point[1]=(3.0+cos(v))*sin(u);
		//		temp_point[2]=sin(v);

		//		TriMesh::Normal d_u,d_v;
		//		d_u[0]=-(3.0+cos(v))*sin(u);
		//		d_u[1]=(3.0+cos(v))*cos(u);
		//		d_u[2]=0.0;
		//		d_v[0]=-cos(u)*sin(v);
		//		d_v[1]=-sin(u)*sin(v);
		//		d_v[2]=cos(v);
		//		TriMesh::Normal temp_normal;
		//		CROSSPROD3(temp_normal, d_u,d_v);
		//		GeometryTools::normalize(temp_normal);
		//		vh_the=cubic_bezier->add_vertex(temp_point);
		//		cubic_bezier->set_normal(vh_the,temp_normal);
		//		vh.push_back(vh_the);

		//		v-=(v1+2.0*PI-v2)/level_num;
		//		if(v<0.0)
		//			v+=2.0*PI;
		//	}
		//}
		////generate each triangle line by line
		//for (int i=0; i<level_num;i++){
		//	int up_index=i*(i+1)/2;
		//	int down_index=(i+1)*(i+2)/2;

		//	for(int j=0; j<i+1;j++){
		//		temp_fh.push_back(vh[up_index]); temp_fh.push_back(vh[down_index]); temp_fh.push_back(vh[down_index+1]);
		//		cubic_bezier->add_face(temp_fh);
		//		temp_fh.clear();
		//		up_index++;
		//		down_index++;
		//	}

		//	up_index=i*(i+1)/2;
		//	down_index=(i+1)*(i+2)/2;
		//	for(int j=0; j<i;j++){
		//		temp_fh.push_back(vh[up_index+1]); temp_fh.push_back(vh[up_index]); temp_fh.push_back(vh[down_index+1]); 
		//		cubic_bezier->add_face(temp_fh);
		//		temp_fh.clear();
		//		up_index++;
		//		down_index++;
		//	}
		//}//generate each triangle line by line
		//vh.clear();




		//u1=(int)(2/12)*PI/12.5;  v1=PI/3.0;
		//u2=(int)(24/12)*PI/12.5;  v2=0.0;
		//u3=(int)(28/12)*PI/12.5;  v3=2.0*PI/3.0;
		////generate each vertex
		//for (int i=0; i<level_num+1;i++){
		//	u=u1+i*(u2-u1)/level_num;
		//	v=v1+i*(v3-v1)/level_num;
		//	for(int j=0; j<i+1;j++){
		//		TriMesh::Point temp_point;
		//		temp_point[0]=(3.0+cos(v))*cos(u);
		//		temp_point[1]=(3.0+cos(v))*sin(u);
		//		temp_point[2]=sin(v);

		//		TriMesh::Normal d_u,d_v;
		//		d_u[0]=-(3.0+cos(v))*sin(u);
		//		d_u[1]=(3.0+cos(v))*cos(u);
		//		d_u[2]=0.0;
		//		d_v[0]=-cos(u)*sin(v);
		//		d_v[1]=-sin(u)*sin(v);
		//		d_v[2]=cos(v);
		//		TriMesh::Normal temp_normal;
		//		CROSSPROD3(temp_normal, d_u,d_v);
		//		GeometryTools::normalize(temp_normal);
		//		vh_the=cubic_bezier->add_vertex(temp_point);
		//		cubic_bezier->set_normal(vh_the,temp_normal);
		//		vh.push_back(vh_the);

		//		v-=(v3-v2)/level_num;
		//	}
		//}
		////generate each triangle line by line
		//for (int i=0; i<level_num;i++){
		//	int up_index=i*(i+1)/2;
		//	int down_index=(i+1)*(i+2)/2;

		//	for(int j=0; j<i+1;j++){
		//		temp_fh.push_back(vh[down_index]); temp_fh.push_back(vh[up_index]); temp_fh.push_back(vh[down_index+1]);
		//		cubic_bezier->add_face(temp_fh);
		//		temp_fh.clear();
		//		up_index++;
		//		down_index++;
		//	}

		//	up_index=i*(i+1)/2;
		//	down_index=(i+1)*(i+2)/2;
		//	for(int j=0; j<i;j++){
		//		temp_fh.push_back(vh[up_index]); temp_fh.push_back(vh[up_index+1]); temp_fh.push_back(vh[down_index+1]); 
		//		cubic_bezier->add_face(temp_fh);
		//		temp_fh.clear();
		//		up_index++;
		//		down_index++;
		//	}
		//}//generate each triangle line by line
		//vh.clear();
}

void database::initalize_cubic_bezier(){
	cubic_bezier =new TriMesh;
	cubic_bezier->request_vertex_normals();
	TriMesh::Point the;
	TriMesh::Normal the_norm;
	TriMesh::VertexHandle vh_the;
	std::vector<TriMesh::VertexHandle> vh;
	std::vector<TriMesh::VertexHandle> temp_fh;
	//TriMesh::VertexIter v_it=themesh->vertices_begin(), v_end=themesh->vertices_end();
	//for ( ; v_it!=v_end; ++v_it){
	//	the=themesh->point(v_it.handle());
	//	vh_the=cubic_bezier->add_vertex(the);
	//	cubic_bezier->set_normal(vh_the, themesh->normal(v_it.handle()));
	//}

	double b1=0.0; 	double b2=0.0; 	double b3=0.0;
	double inc=1.0/3.0;
	double coeff[10][3];
	double normCoeff[6][3];

	TriMesh::FaceIter f_it=themesh->faces_begin(), f_end=themesh->faces_end();
	for ( ; f_it!=f_end; ++f_it)	{
		
		//themesh->Cubic_Bezier_Parameters(f_it.handle(), coeff, normCoeff);

		vh.clear();
		b1=1.0; b2=0.0;b3=0.0;
		themesh->para_to_point_PN_triangle(f_it.handle(),  b1, b2,  b3, the, the_norm);
		vh_the=cubic_bezier->add_vertex(the);
		cubic_bezier->set_normal(vh_the,the_norm);
		vh.push_back(vh_the);
		///////////////////////////////////////////////////////////////////////////////////
		b1=0.0; b2=1.0;
		themesh->para_to_point_PN_triangle(f_it.handle(),  b1, b2,  b3, the, the_norm);
		vh_the=cubic_bezier->add_vertex(the);
		cubic_bezier->set_normal(vh_the,the_norm);
		vh.push_back(vh_the);
		///////////////////////////////////////////////////////////////////////////////////
		b2=0.0;b3=1.0;
		themesh->para_to_point_PN_triangle(f_it.handle(),  b1, b2,  b3, the, the_norm);
		vh_the=cubic_bezier->add_vertex(the);
		cubic_bezier->set_normal(vh_the,the_norm);
		vh.push_back(vh_the);

		b1=0.0; 		b2=1.0-inc; 		b3=inc;
		for (int j=1; j<4; j++) { //For each row of vertices
			b1 = inc*j;
			b3 = 1.0f - b1 - b2;	
			for (int k=0; k<j+1; k++) //For each vertex in that row
			{
				if( (1.0-b1)<0.00000000001||(1.0-b3)<0.00000000001)
				{
					b1 -= inc; 	 b3 += inc;		continue;
				}
				themesh->para_to_point_PN_triangle(f_it.handle(),  b1, b2,  b3, the, the_norm);
				//themesh->EvaluateBezierPoint(coeff, normCoeff, b1,b2,b3, the, the_norm);
				vh_the=cubic_bezier->add_vertex(the);
				cubic_bezier->set_normal(vh_the,the_norm);
				vh.push_back(vh_the);				
				b1 -= inc;
				b3 += inc;
			}			
			b2 -= inc;
		}

		//row 1
		temp_fh.clear();
		temp_fh.push_back(vh[3]); temp_fh.push_back(vh[1]); temp_fh.push_back(vh[4]);
		cubic_bezier->add_face(temp_fh);
		temp_fh.clear();
		//row 2
		temp_fh.push_back(vh[5]); temp_fh.push_back(vh[3]); temp_fh.push_back(vh[6]);
		cubic_bezier->add_face(temp_fh);
		temp_fh.clear();
		temp_fh.push_back(vh[6]); temp_fh.push_back(vh[3]); temp_fh.push_back(vh[4]);
		cubic_bezier->add_face(temp_fh);
		temp_fh.clear();
		temp_fh.push_back(vh[6]); temp_fh.push_back(vh[4]); temp_fh.push_back(vh[7]);
		cubic_bezier->add_face(temp_fh);
		temp_fh.clear();
		//row 3
		temp_fh.push_back(vh[0]); temp_fh.push_back(vh[5]); temp_fh.push_back(vh[8]);
		cubic_bezier->add_face(temp_fh);
		temp_fh.clear();
		temp_fh.push_back(vh[8]); temp_fh.push_back(vh[5]); temp_fh.push_back(vh[6]);
		cubic_bezier->add_face(temp_fh);
		temp_fh.clear();
		temp_fh.push_back(vh[8]); temp_fh.push_back(vh[6]); temp_fh.push_back(vh[9]);
		cubic_bezier->add_face(temp_fh);
		temp_fh.clear();
		temp_fh.push_back(vh[9]); temp_fh.push_back(vh[6]); temp_fh.push_back(vh[7]);
		cubic_bezier->add_face(temp_fh);
		temp_fh.clear();
		temp_fh.push_back(vh[9]); temp_fh.push_back(vh[7]); temp_fh.push_back(vh[2]);
		cubic_bezier->add_face(temp_fh);
		temp_fh.clear();

	}
	vh.clear();
}

void database::clear(){
	if(textureimage){
		delete [] textureimage;
		textureimage=NULL;
	}
	if(themesh){
		delete [] themesh;
		themesh=NULL;
	}

	if(cubic_bezier){
		delete [] cubic_bezier;
		cubic_bezier=NULL;
	}

}

void database::save(){
	
}

void database::read(){
	
}
