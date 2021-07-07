from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

import numpy as np
from collections import namedtuple
import os

#E: [MPa], nu: [-], rho: [Mg/mm^3]
YarnType = namedtuple('YarnType',['E','nu','rho_density'])
yarn_prop_paper = YarnType(1400.0,0.2,1.306e-09)
yarn_prop_pink = YarnType(51.5283,0.2,8.5532e-10)
#todo: recalc yarn_prop_pink E and rho w/ lower r

r_fake = 0.08
#r_pink = 0.8059 #this is what was measured vs below is what worked
#r_pink = 0.255
r_pink = 0.245

StitchSize = namedtuple('StitchSize',['lam', 'w', 'gamma', 'CO', 'h', 'delta'])
stitch_size_paper = StitchSize(3.73,1.15,0.24,0.71,1.14,0.41)
stitch_size_5 = StitchSize(11.7123,3.0417,0.6045,2.5345,4.2130,2.0000)

def run_inp(jname, num_threads = 4):
	if os.path.exists(jname + '.inp'):
		# j = mdb.JobFromInputFile(atTime=None, explicitPrecision=DOUBLE_PLUS_PACK, 
		# 	getMemoryFromAnalysis=True, inputFileName=jname + '.inp'
		# 	, memory=90, memoryUnits=PERCENTAGE, multiprocessingMode=DEFAULT, name=jname, 
		# 	nodalOutputPrecision=SINGLE, numCpus=num_threads,numDomains = num_threads numGPUs=0, 
		# 	queue=None, resultsFormat=ODB, scratch='', type=ANALYSIS, userSubroutine=''
		# 	, waitHours=0, waitMinutes=0)

		#a = 'C:\\Users\\read8\\Google Drive\\research\\abaqus modeling\\job-test.inp'
		j = mdb.JobFromInputFile(activateLoadBalancing=False, atTime=None, 
			explicitPrecision=DOUBLE_PLUS_PACK, inputFileName=jname, 
			memory=90, memoryUnits=PERCENTAGE, multiprocessingMode=DEFAULT, name=
			jname, nodalOutputPrecision=SINGLE, numCpus=num_threads, numDomains=num_threads, 
			parallelizationMethodExplicit=DOMAIN, queue=None, resultsFormat=ODB, 
			scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

		j.submit(consistencyChecking = OFF)
		j.waitForCompletion()
	#end
#end

def get_fd(jname):
	#post processing
	file_name = jname + '.odb'
	path_j = file_name
	odb = openOdb(path = path_j)
	step = odb.steps['Step-1']
	his_name = step.historyRegions.keys()[1]
	his_region = step.historyRegions[his_name]

	out_his_RF2 = his_region.historyOutputs['RF2'].data
	out_his_u2 = his_region.historyOutputs['U2'].data
	odb.close()

	rf2 = [out_his_RF2[idx][1] for idx in range(len(out_his_RF2))]
	u2 = [out_his_u2[idx][1] for idx in range(len(out_his_u2))]

	f = open('_fd_' + str(jname) + '.txt','w')
	for idx in range(len(rf2)):
	    f.write(str(u2[idx]) + ' ' + str(rf2[idx]) + '\n')
	#end
	f.close()
#end

def delete_extra_files(jname):
	#delete extra files
	extra_file_extensions = ['.com', '.dat', '.inp', '.log','.msg','.prt','.sim','.sta']
	for i in range(len(extra_file_extensions)):
		if os.path.exists(jname + extra_file_extensions[i]):
			os.remove(jname + extra_file_extensions[i])
		#end
	#end
#end

def printAB(string_):
	print >> sys.__stdout__, string_
#end

class knit_beam:
	cut_point = 1.0
	cut_size = 0.5

	#set value
	fric_coeff = 0.3

	def __init__(self, stitchProp,matType,r):
		self.lam = stitchProp.lam
		self.w = stitchProp.w
		self.gamma = stitchProp.gamma
		self.CO = stitchProp.CO
		self.h = stitchProp.h
		self.delta = stitchProp.delta

		self.E = matType.E
		self.nu = matType.nu
		self.rho_density = matType.rho_density

		self.r = r
	#end
	def xyz_swept(self,s):
		lam = self.lam; gamma = self.gamma; w = self.w;
		CO = self.CO; h = self.h; delta = self.delta;

		x = s/lam * w + gamma*np.sin(4*np.pi*s/lam)
		y = CO + h/2 * np.cos(2*np.pi*s/lam)
		z = delta/2 * np.cos(4*np.pi*s/lam)

		return (x,y,z)
	#end

	def change_friction(self, fric_coeff_):
		self.fric_coeff = fric_coeff_
	#end

	def make_cut_knit(self, num_cycles, num_rows, cae_name, cut_point_, cut_size_):
		self.cut_size = cut_size_
		self.cut_point = cut_point_
		cut_row = np.ceil(num_rows/2.0) #one indexed (sry!!)
		self.make_knit(num_cycles, num_rows, cae_name, cut_row)
	#end

	def make_perf_knit(self, num_cycles, num_rows, cae_name):
		self.make_knit(num_cycles, num_rows, cae_name)
	#end

	def make_knit(self,num_cycles, num_rows, cae_name, cut_row = -1):
		Mdb()
		m = mdb.models['Model-1']
		p = m.Part(name='Part-knit', dimensionality=THREE_D,type=DEFORMABLE_BODY)

		y_offset = self.CO; h = self.h; r = self.r; delta = self.delta; lam = self.lam

		E = self.E; nu = self.nu; rho_density = self.rho_density; fric_coeff = self.fric_coeff;

		cut_size = self.cut_size
		cut_point = self.cut_point

		tol_bb = 2*r

		num_cycles = float(num_cycles)
		num_rows = int(num_rows)

		disp_top = 11.0 #mm

		#when num_cycles was 5.0 it was 1.5,2.5
		t_cut_1 = cut_point - float(cut_size/2)
		t_cut_2 = cut_point + float(cut_size/2)

		t_len = 2000
		t_len_cut = 1000


		xyz_cut_left = []
		xyz_cut_right = []
		t_all = np.linspace(0,(num_cycles - 1)*lam,num=t_len)
		t_cut_left = np.linspace(0,t_cut_1*lam, num=t_len_cut)
		t_cut_right = np.linspace(t_cut_2*lam, (num_cycles - 1)*lam, num=t_len_cut)

		v0 = self.xyz_swept(0)
		vf = self.xyz_swept(t_all[t_len - 1])

		xyz = [None] * t_len
		xyz_left = [None] * t_len_cut
		xyz_right = [None] * t_len_cut

		#geometry/sketching
		for i in range(num_rows):
			if i == (cut_row - 1):
				for j in range(t_len_cut):
					pt_left = self.xyz_swept(t_cut_left[j])
					pt_right = self.xyz_swept(t_cut_right[j])

					y_L = pt_left[1] + i*y_offset
					y_R = pt_right[1] + i*y_offset

					xyz_left[j] = (pt_left[0], y_L, pt_left[2])
					xyz_right[j] = (pt_right[0], y_R, pt_right[2])
				#end
				p.WireSpline(points=xyz_left, meshable=ON,smoothClosedSpline=ON)
				p.WireSpline(points=xyz_right, meshable=ON,smoothClosedSpline=ON)
			#end
			else:
				for j in range(len(t_all)):
					point_cur = self.xyz_swept(t_all[j])
					x_cur = point_cur[0]
					y_cur = point_cur[1]
					z_cur = point_cur[2]

					point_et = (x_cur, y_cur + i*y_offset, z_cur)
					xyz[j] = point_et
				#end
				p.WireSpline(points=xyz, meshable=ON,smoothClosedSpline=ON)
			#end
		#end

		#assembly
		aa = m.rootAssembly
		i_all = aa.Instance(dependent=ON, name='Part-knit', part=p)

		#material
		#elastic table has units of MPa
		#density table has units of Mg/mm^3
		mat_yarn = m.Material(name='mat-yarn')
		# mat_yarn.Elastic(table=((74900.0, 1400.0, 
		#     1400.0, 0.1, 0.1, 0.35, 440.0, 440.0, 440.0), ), type=
		#     ENGINEERING_CONSTANTS)
		mat_yarn.Elastic(table=((E, nu), ))
		mat_yarn.Density(table=((rho_density, ), ))

		#set knit
		edges_yarn = p.edges.getByBoundingBox(xMin = 0-tol_bb, 
		    xMax = vf[0]+tol_bb,
		    yMin = y_offset - h/2 - tol_bb, 
		    yMax = num_rows*y_offset + h/2 + tol_bb,
		    zMin = -delta/2 - tol_bb, 
		    zMax = delta/2 + tol_bb,)
		set_yarn = p.Set(edges = edges_yarn, name='Set-yarn')

		#section
		m.CircularProfile(name='Profile-yarn', r=r)
		m.BeamSection(consistentMassMatrix=False, integration=DURING_ANALYSIS, material='mat-yarn', 
			name='Section-yarn', poissonRatio=0.0, profile='Profile-yarn', temperatureVar=LINEAR)
		p.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=set_yarn, 
			sectionName='Section-yarn', thicknessAssignment=FROM_SECTION)
		p.assignBeamSectionOrientation(method=N1_COSINES, n1=(0.0, 0.0, -1.0), region=set_yarn)

		#mesh
		p.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size = lam/16.0)
		p.generateMesh()
		aa.regenerate()

		#set: roller bd
		vert_left = i_all.vertices.getByBoundingBox(xMin = 0 - tol_bb, xMax = tol_bb,
			yMin = y_offset - h/2 - tol_bb, yMax = num_rows*y_offset + h/2 + tol_bb,
			zMin = -delta/2 - tol_bb, zMax = delta/2 + tol_bb)

		vert_right = i_all.vertices.getByBoundingBox(xMin = vf[0] - tol_bb, xMax = vf[0] + tol_bb,
			yMin = y_offset - h/2 - tol_bb, yMax = num_rows*y_offset + h/2 + tol_bb,
			zMin = -delta/2 - tol_bb, zMax = delta/2 + tol_bb)

		set_roller_edge = aa.Set(name='Set-roller-bd', vertices=[vert_left,vert_right])

		#set: top/bottom
		y_con = 0.5*h

		nodes_top = i_all.nodes.getByBoundingBox(xMin = 0-tol_bb, xMax = vf[0]+tol_bb,
		    yMin = (num_rows*y_offset + r + h/2) - y_con, yMax = num_rows*y_offset + r + h/2 + tol_bb,
		    zMin =  -r -delta/2 -tol_bb, zMax = r + delta/2 +tol_bb,)
		set_top = aa.Set(name = 'Set-top-bc', nodes = nodes_top)

		nodes_bottom = i_all.nodes.getByBoundingBox(xMin = 0-tol_bb, xMax = vf[0]+tol_bb,
		    yMin = y_offset - r - h/2 - tol_bb, yMax = (y_offset - r - h/2) + y_con,
		    zMin = -r -delta/2 - tol_bb, zMax = r + delta/2 + tol_bb,)
		set_bottom = aa.Set(name = 'Set-bottom-bc', nodes = nodes_bottom)

		#step
		m.ExplicitDynamicsStep(improvedDtMethod=ON, name='Step-1', previous='Initial')

		#bd cond: roller/u1 anti?symmetric
		# zero_face = (0,(i+1)*y_offset + h/2,z_offset + delta/2)
	 #    face_list.append(i_all.faces.findAt((zero_face,)))
	 #    end_face = i_all.faces.getByBoundingSphere(center = (vf[0],(i+1)*y_offset + h/2,z_offset + delta/2), radius=tol)
	 # end_face = i_all.faces.getByBoundingSphere(center = (vf[0],(i+1)*y_offset + h/2,z_offset + delta/2), radius=tol)
		for i in range(num_rows):
			vert_cur_left = i_all.vertices.getByBoundingSphere(center = (v0[0],(i+1)*y_offset + h/2, delta/2), 
				radius = tol_bb)
			vert_cur_right = i_all.vertices.getByBoundingSphere(center = (vf[0],(i+1)*y_offset + h/2, delta/2), 
				radius = tol_bb)

			set_nameL = 'Set-nodeL-' + str(i+1)
			set_nameR = 'Set-nodeR-' + str(i+1)
			vert_left_set = aa.Set(name = set_nameL, vertices = vert_cur_left)
			vert_right_set = aa.Set(name = set_nameR, vertices = vert_cur_right)

			eqn_name_u1 = 'u1-row-' + str(i+1)
			eqn_name_u2 = 'u2-row-' + str(i+1)
			eqn_name_u3 = 'u3-row-' + str(i+1)
			eqn_name_ur1 = 'ur1-row-' + str(i+1)
			eqn_name_ur2 = 'ur2-row-' + str(i+1)
			eqn_name_ur3 = 'ur3-row-' + str(i+1)

			m.Equation(name = eqn_name_u1, terms = ((1.0, set_nameR, 1), (1.0, set_nameL, 1)))
			m.Equation(name = eqn_name_u2, terms = ((1.0, set_nameR, 2), (-1.0, set_nameL, 2)))
			m.Equation(name = eqn_name_u3, terms = ((1.0, set_nameR, 3), (-1.0, set_nameL, 3)))
			# m.Equation(name = eqn_name_ur1, terms = ((1.0, set_nameR, 3), (-1.0, set_nameL, 4)))
			# m.Equation(name = eqn_name_ur2, terms = ((1.0, set_nameR, 3), (-1.0, set_nameL, 5)))
			# m.Equation(name = eqn_name_ur3, terms = ((1.0, set_nameR, 3), (-1.0, set_nameL, 6)))
			#m.Equation(name='Constraint-top', terms=((1.0, 'Set-top-bc', 2), (-1.0, 'Set-RP', 2)))
		#end
		
		# m.DisplacementBC(amplitude=UNSET, createStepName='Initial', distributionType=UNIFORM, fieldName='', 
		# 	localCsys=None, name='BC-roller-edge', region=set_roller_edge, u1=UNSET, u2=UNSET,
		# 	 u3=UNSET, ur1=SET, ur2=SET, ur3=SET)
		#what we actually want is that each pair on L/R has periodic BC
		#u1_L + u1_R = 0 they should be opposite??????
		#u2_L - u2_R = 0 they should be the same
		#u3_L - u3_R = 0
		#ur1 is same, ur2 is same [??], ur3 is opp

		#what if I do something kinda simple and restrict u1 to an eqn and let everthing else just be set/unset as before

		#bd cond:top/bottom:
		ref_pt = aa.ReferencePoint(point=(1.0, 1.0, 1.0))
		set_rp = aa.Set(name='Set-RP', referencePoints=(aa.referencePoints[ref_pt.id], ))
		m.Equation(name='Constraint-top', terms=((1.0, 'Set-top-bc', 2), (-1.0, 'Set-RP', 2)))

		m.EncastreBC(createStepName='Initial', localCsys=None, name='BC-bottom', region=set_bottom)
		m.EncastreBC(createStepName='Initial', localCsys=None, name='BC-top-initial', region=set_top)

		m.TabularAmplitude(data=((0.0, 0.0), (1.0, 1.0)), name='Amp-ramp', smooth=SOLVER_DEFAULT, timeSpan=STEP)
		m.DisplacementBC(amplitude='Amp-ramp', createStepName='Step-1'
		    , distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='BC-top', region=set_rp, 
		    u1=UNSET, u2=disp_top, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

		#history output
		m.HistoryOutputRequest(createStepName='Step-1', name='H-Output-rp', rebar=EXCLUDE, region=
		    aa.sets['Set-RP'], sectionPoints=DEFAULT, variables=('U2', 'RF2'))

		#contact prop (total guesswork here)
		m.ContactProperty('IntProp-contact')
		m.interactionProperties['IntProp-contact'].TangentialBehavior(
		    dependencies=0, directionality=ISOTROPIC, elasticSlipStiffness=None, 
		    formulation=PENALTY, fraction=0.005, maximumElasticSlip=FRACTION, 
		    pressureDependency=OFF, shearStressLimit=None, slipRateDependency=OFF, 
		    table=((fric_coeff, ), ), temperatureDependency=OFF)
		m.interactionProperties['IntProp-contact'].NormalBehavior(
		    allowSeparation=ON, constraintEnforcementMethod=DEFAULT, 
		    pressureOverclosure=HARD)

		#contact :(((((
		int_contact = m.ContactExp(createStepName='Initial', name='Int-contact')
		int_contact.includedPairs.setValuesInStep(stepName='Initial', useAllstar=ON)
		int_contact.contactPropertyAssignments.appendInStep(
		    assignments=((GLOBAL, SELF, 'IntProp-contact'), ), stepName='Initial')

		#field output request- nicer movies
		m.fieldOutputRequests['F-Output-1'].setValues(numIntervals=72)

		mdb.saveAs(cae_name + '.cae')
	#end

	def make_inp(self, jname, disp_top, num_threads = 4):
		m = mdb.models['Model-1']
		disp_top = float(disp_top)
		m.boundaryConditions['BC-top'].setValues(u2 = disp_top)

		j = mdb.Job(activateLoadBalancing=False, atTime=None, contactPrint=OFF, 
			description='', echoPrint=OFF, explicitPrecision=DOUBLE_PLUS_PACK, 
			historyPrint=OFF, memory=90, memoryUnits=PERCENTAGE, model='Model-1', 
			modelPrint=OFF, multiprocessingMode=DEFAULT, name = jname, 
			nodalOutputPrecision=SINGLE, numCpus=num_threads, numDomains=num_threads, 
			parallelizationMethodExplicit=DOMAIN, queue=None, resultsFormat=ODB, 
			scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

		j.writeInput()

	#end

#end;

class knit_3d:
	cut_point = 1.0
	cut_size = 0.5

	#set value
	fric_coeff = 0.3

	def __init__(self, stitchProp,matType,r):
		self.lam = stitchProp.lam
		self.w = stitchProp.w
		self.gamma = stitchProp.gamma
		self.CO = stitchProp.CO
		self.h = stitchProp.h
		self.delta = stitchProp.delta

		self.E = matType.E
		self.nu = matType.nu
		self.rho_density = matType.rho_density

		self.r = r
	#end
	def xyz_swept(self,s):
		lam = self.lam; gamma = self.gamma; w = self.w;
		CO = self.CO; h = self.h; delta = self.delta;

		x = s/lam * w + gamma*np.sin(4*np.pi*s/lam)
		y = CO + h/2 * np.cos(2*np.pi*s/lam)
		z = delta/2 * np.cos(4*np.pi*s/lam)

		return (x,y,z)
	#end

	def make_cut_knit(self, num_cycles, num_rows, cae_name, cut_point_, cut_size_):
		self.cut_size = float(cut_size_)
		self.cut_point = float(cut_point_)
		cut_row = np.ceil(num_rows/2.0) #one indexed (sry!!)
		self.make_knit(num_cycles, num_rows, cae_name, cut_row)
	#end

	def make_perf_knit(self, num_cycles, num_rows, cae_name):
		self.make_knit(num_cycles, num_rows, cae_name)
	#end

	def make_knit(self, num_cycles, num_rows, cae_name, cut_row = -1):
		Mdb()
		m = mdb.models['Model-1']
		p = m.Part(name='Part-knit', dimensionality=THREE_D,type=DEFORMABLE_BODY)

		y_offset = self.CO; h = self.h; r = self.r; delta = self.delta; lam = self.lam

		E = self.E; nu = self.nu; rho_density = self.rho_density; fric_coeff = self.fric_coeff;

		cut_size = self.cut_size
		cut_point = self.cut_point

		#BCs top/bottom
		disp_top = 11.0

		tol_bb = 2*r

		num_cycles = float(num_cycles)
		num_rows = int(num_rows)

		z_offset = 2*delta

		t_cut_1 = cut_point - float(cut_size/2)
		t_cut_2 = cut_point + float(cut_size/2)

		t_len = 2000

		norm_len = (num_cycles - 1)
		left_len = t_cut_1
		right_len = (num_cycles - 1) - t_cut_2

		t_len_left = int(t_len * left_len/norm_len)
		t_len_right = int(t_len * right_len/norm_len)


		#t_len_cut = 1000


		xyz=[]
		xyz_cut_left = []
		xyz_cut_right = []
		t_all = np.linspace(0,(num_cycles - 1)*lam,num=t_len)
		t_cut_left = np.linspace(0,t_cut_1*lam, num=t_len_left)
		t_cut_right = np.linspace(t_cut_2*lam, (num_cycles - 1)*lam, num=t_len_right)

		v0 = self.xyz_swept(0)
		vf = self.xyz_swept(t_all[t_len - 1])

		v0_L = (v0[0],v0[1] + 2*y_offset, v0[2])
		vf_R = (vf[0],vf[1] + 3*y_offset, vf[2])

		for i in range(len(t_all)):
		    xyz.append(self.xyz_swept(t_all[i]))
		#end      

		p.WireSpline(points=xyz, meshable=OFF,smoothClosedSpline=ON)

		for i in range(len(t_cut_left)):
			left_pt = self.xyz_swept(t_cut_left[i])
			right_pt = self.xyz_swept(t_cut_right[i])

			left_pt = (left_pt[0], left_pt[1]+2*y_offset, left_pt[2])
			right_pt = (right_pt[0], right_pt[1] + 3*y_offset, right_pt[2])

			xyz_cut_left.append(left_pt)
			xyz_cut_right.append(right_pt)
		#end


		p.WireSpline(points=xyz_cut_left, meshable=OFF,smoothClosedSpline=ON)
		p.WireSpline(points=xyz_cut_right, meshable=OFF,smoothClosedSpline=ON)


		#datums for p
		p.DatumPointByCoordinate(coords=(0.0, 0.0, 0.0))
		yz_plane = p.DatumPlaneByPrincipalPlane(offset=0.0, principalPlane=YZPLANE)
		yz_plane_offset = p.DatumPlaneByPrincipalPlane(offset=vf[0], principalPlane=YZPLANE)

		y_axis = p.DatumAxisByPrincipalAxis(principalAxis=YAXIS)


		size = 8*h

		for i in range(num_rows):
			if i == (cut_row - 1):
				#left side:
				sweep_along_L = p.edges.findAt((v0_L, ))
				set_L = p.Set(edges = sweep_along_L, name='Set-L-edge')
				s_row_cur = m.ConstrainedSketch(gridSpacing=0.21, name='s_row_cur', 
					sheetSize=size, transform=
					p.MakeSketchTransform(sketchPlane=p.datums[yz_plane.id], sketchPlaneSide=SIDE1, 
					sketchUpEdge=p.datums[y_axis.id], sketchOrientation=BOTTOM, origin=(0.0, 0.0, z_offset)))
				s_row_cur.ConstructionLine(point1=(-size/2, 0.0), point2=(size/2, 0.0))
				s_row_cur.ConstructionLine(point1=(0.0, -size/2), point2=(0.0, size/2))

				center_circle = (v0[1] + i*y_offset,v0[2])
				perim_circle = (v0[1] + i*y_offset + r, v0[2])
				s_row_cur.CircleByCenterPerimeter(center=center_circle, point1=perim_circle)
				p.SolidSweep(path=sweep_along_L, profile=s_row_cur, 
					sketchOrientation=BOTTOM, sketchPlane=p.datums[yz_plane.id], sketchUpEdge=p.datums[y_axis.id])
				del s_row_cur

				#right side:
				sweep_along_R = p.edges.findAt((vf_R, ))
				set_R = p.Set(edges = sweep_along_R, name='Set-R-edge')
				s_row_cur = m.ConstrainedSketch(gridSpacing=0.21, name='s_row_cur', 
					sheetSize=size, transform=
					p.MakeSketchTransform(sketchPlane=p.datums[yz_plane_offset.id], sketchPlaneSide=SIDE1, 
					sketchUpEdge=p.datums[y_axis.id], sketchOrientation=BOTTOM, origin=(0.0, 0.0, z_offset)))
				s_row_cur.ConstructionLine(point1=(-size/2, 0.0), point2=(size/2, 0.0))
				s_row_cur.ConstructionLine(point1=(0.0, -size/2), point2=(0.0, size/2))

				center_circle = (v0[1] + i*y_offset,v0[2])
				perim_circle = (v0[1] + i*y_offset + r, v0[2])
				s_row_cur.CircleByCenterPerimeter(center=center_circle, point1=perim_circle)
				p.SolidSweep(flipSweepDirection=ON, path=sweep_along_R, profile=s_row_cur, 
					sketchOrientation=BOTTOM, sketchPlane=p.datums[yz_plane_offset.id], sketchUpEdge=p.datums[y_axis.id])
				del s_row_cur
			else:
				sweep_along = p.edges.findAt((v0, ))
				s_row_cur = m.ConstrainedSketch(gridSpacing=0.21, name='s_row_cur', 
					sheetSize=size, transform=
					p.MakeSketchTransform(sketchPlane=p.datums[yz_plane.id], sketchPlaneSide=SIDE1, 
					sketchUpEdge=p.datums[y_axis.id], sketchOrientation=BOTTOM, origin=(0.0, 0.0, z_offset)))
				s_row_cur.ConstructionLine(point1=(-size/2, 0.0), point2=(size/2, 0.0))
				s_row_cur.ConstructionLine(point1=(0.0, -size/2), point2=(0.0, size/2))
				#p.projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=s_row_cur)
				center_circle = (v0[1] + i*y_offset,v0[2])
				perim_circle = (v0[1] + i*y_offset + r, v0[2])
				s_row_cur.CircleByCenterPerimeter(center=center_circle, point1=perim_circle)
				p.SolidSweep(path=sweep_along, profile=s_row_cur, 
					sketchOrientation=BOTTOM, sketchPlane=p.datums[yz_plane.id], sketchUpEdge=p.datums[y_axis.id])
				del s_row_cur
			#end
		#end


		#assembly
		aa = m.rootAssembly
		aa.DatumCsysByDefault(CARTESIAN)
		i_all = aa.Instance(dependent=ON, name='Part-knit-1', part=p)


		#set roller bd: edge
		#tol = 0.1
		face_list = []
		for i in range(num_rows):
		    idx = 2*i;
		    zero_face = (0,(i+1)*y_offset + h/2,z_offset + delta/2)
		    face_list.append(i_all.faces.findAt((zero_face,)))
		    end_face = i_all.faces.getByBoundingSphere(center = (vf[0],(i+1)*y_offset + h/2,z_offset + delta/2), radius=tol_bb)
		    face_list.append(end_face)

		    name_left_cur = 'Set-faceL-' + str(int(i+1))
		    name_right_cur = 'Set-faceR-' + str(int(i+1))
		    name_ref_pt_cur = 'Set-RP-bd' + str(int(i+1))
		    ref_pt_bd_cur = aa.ReferencePoint(point=(2.0, i*y_offset, z_offset))

		    eqn_nameL = 'u1-rowL-' + str(int(i+1))
		    eqn_nameR = 'u1-rowR-' + str(int(i+1))
		    face_left = aa.Set(faces = i_all.faces.findAt((zero_face,)), name = name_left_cur)
		    face_right = aa.Set(faces = end_face, name = name_right_cur)
		    set_ref_pt_bd_cur = aa.Set(name = name_ref_pt_cur, referencePoints = (aa.referencePoints[ref_pt_bd_cur.id],))

		    m.Equation(name = eqn_nameL, terms = ((1.0, name_left_cur, 1), (1.0, name_ref_pt_cur, 1)))
		    m.Equation(name = eqn_nameR, terms = ((1.0, name_right_cur, 1), (-1.0, name_ref_pt_cur, 1)))
		#end
		set_roller_edge = aa.Set(faces = face_list, name = 'Set-roller-bd')


		#bd roller
		m.DisplacementBC(amplitude=UNSET, createStepName='Initial', 
		    distributionType=UNIFORM, fieldName='', localCsys=None, name=
		    'BC-roller-edge', region=
		    set_roller_edge, u1=UNSET, u2=UNSET, 
		    u3=SET, ur1=SET, ur2=SET, ur3=SET)

		#set for wire (booooo)
		edges_wire= p.edges.getByBoundingBox(xMin = 0-tol_bb, xMax = vf[0]+tol_bb,yMin = y_offset - h/2 -tol_bb, 
		    yMax = 4*y_offset + h/2 +tol_bb,zMin = -delta/2 -tol_bb, zMax = delta/2 +tol_bb,)
		set_wire = p.Set(edges = edges_wire, name='Set-wire')
		m.EncastreBC(createStepName='Initial', localCsys=None, 
		    name='BC-wire', region=i_all.sets['Set-wire'])

		#step
		m.ExplicitDynamicsStep(improvedDtMethod=ON, name='Step-1', previous='Initial')

		#material
		#elastic table has units of MPa
		#density table has units of Mg/mm^3
		mat_yarn = m.Material(name='mat-yarn')
		# mat_yarn.Elastic(table=((74900.0, 1400.0, 
		#     1400.0, 0.1, 0.1, 0.35, 440.0, 440.0, 440.0), ), type=
		#     ENGINEERING_CONSTANTS)
		mat_yarn.Elastic(table=((1400, 0.2), ))
		mat_yarn.Density(table=((1.306e-09, ), )) #oops should be 1e-9!!! todo for future work

		#material fake
		m.Material(name='Material-fake')
		m.materials['Material-fake'].Elastic(table=((1.0, 0.0), ))
		m.materials['Material-fake'].Density(table=((1e-09, ), ))

		#section real
		m.HomogeneousSolidSection(material='mat-yarn', name='Section-yarn', thickness=None)

		cells_yarn = p.cells.getByBoundingBox(xMin = 0-tol_bb, 
		    xMax = vf[0]+tol_bb,
		    yMin = y_offset - r - h/2 -tol_bb, 
		    yMax = num_rows*y_offset + r + h/2 +tol_bb,
		    zMin = z_offset - r -delta/2 -tol_bb, 
		    zMax = z_offset + r + delta/2 +tol_bb,)
		set_yarn = p.Set(cells = cells_yarn, name='Set-yarn')


		p.SectionAssignment(offset=0.0, 
		    offsetField='', offsetType=MIDDLE_SURFACE, region=set_yarn, 
		    sectionName='Section-yarn', thicknessAssignment=FROM_SECTION)

		#section fake
		m.TrussSection(area=1.0, material='Material-fake', name='Section-wire')
		p.SectionAssignment(offset=0.0, 
		    offsetField='', offsetType=MIDDLE_SURFACE, region=
		    p.sets['Set-wire'], sectionName=
		    'Section-wire', thicknessAssignment=FROM_SECTION)

		#contact prop (total guesswork here)
		m.ContactProperty('IntProp-contact')
		m.interactionProperties['IntProp-contact'].TangentialBehavior(
		    dependencies=0, directionality=ISOTROPIC, elasticSlipStiffness=None, 
		    formulation=PENALTY, fraction=0.005, maximumElasticSlip=FRACTION, 
		    pressureDependency=OFF, shearStressLimit=None, slipRateDependency=OFF, 
		    table=((0.3, ), ), temperatureDependency=OFF)
		m.interactionProperties['IntProp-contact'].NormalBehavior(
		    allowSeparation=ON, constraintEnforcementMethod=DEFAULT, 
		    pressureOverclosure=HARD)

		#contact
		m.ContactExp(createStepName='Initial', name='Int-contact')
		m.interactions['Int-contact'].includedPairs.setValuesInStep(
		    stepName='Initial', useAllstar=ON)
		m.interactions['Int-contact'].contactPropertyAssignments.appendInStep(
		    assignments=((GLOBAL, SELF, 'IntProp-contact'), ), stepName='Initial')

		#mesh
		p.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=lam/16)
		p.setElementType(elemTypes=(ElemType(
		    elemCode=T3D2, elemLibrary=STANDARD), ), regions=set_wire)
		p.setElementType(elemTypes=(ElemType(
		    elemCode=C3D8R, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
		    kinematicSplit=AVERAGE_STRAIN, hourglassControl=DEFAULT, 
		    distortionControl=DEFAULT), ElemType(elemCode=C3D6, elemLibrary=STANDARD), 
		    ElemType(elemCode=C3D4, elemLibrary=STANDARD)), regions=set_yarn)
		p.generateMesh()
		aa.regenerate()


		#sets/bc top/bottom
		y_con = 0.5

		nodes_top = i_all.nodes.getByBoundingBox(xMin = 0-tol_bb, 
		    xMax = vf[0]+tol_bb,
		    yMin = (num_rows*y_offset + r + h/2) - y_con, 
		    yMax = num_rows*y_offset + r + h/2 + tol_bb,
		    zMin = z_offset - r -delta/2 -tol_bb, 
		    zMax = z_offset + r + delta/2 +tol_bb,)

		aa.Set(name = 'Set-top-bc', nodes = nodes_top)

		nodes_bottom = i_all.nodes.getByBoundingBox(xMin = 0-tol_bb, 
		    xMax = vf[0]+tol_bb,
		    yMin = y_offset - r - h/2 -tol_bb, 
		    yMax = (y_offset - r - h/2) + y_con,
		    zMin = z_offset - r -delta/2 -tol_bb, 
		    zMax = z_offset + r + delta/2 +tol_bb,)

		aa.Set(name = 'Set-bottom-bc', nodes = nodes_bottom)


		#create reference point
		ref_pt = aa.ReferencePoint(point=(1.0, 1.0, 1.0))
		aa.Set(name='Set-RP', referencePoints=(aa.referencePoints[ref_pt.id], ))
		m.Equation(name='Constraint-top', terms=((1.0, 'Set-top-bc', 2), (-1.0, 'Set-RP', 2)))

		#history output
		m.HistoryOutputRequest(createStepName='Step-1', name='H-Output-rp', rebar=EXCLUDE, region=
		    aa.sets['Set-RP'], sectionPoints=DEFAULT, variables=('U2', 'RF2'))

		m.DisplacementBC(amplitude=UNSET, createStepName='Initial', 
		    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-bottom', 
		    region=aa.sets['Set-bottom-bc'], u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET)

		m.TabularAmplitude(data=((0.0, 0.0), (1.0, 1.0)), name='Amp-ramp', smooth=SOLVER_DEFAULT, timeSpan=STEP)
		m.DisplacementBC(amplitude='Amp-ramp', createStepName='Step-1'
		    , distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='BC-top', region=aa.sets['Set-RP'], 
		    u1=0, u2=disp_top, u3=0, ur1=0, ur2=0, ur3=0)

		m.fieldOutputRequests['F-Output-1'].setValues(numIntervals=72)

		mdb.saveAs(cae_name + '.cae')
	#end

	def make_inp(self, jname, disp_top, num_threads = 4):
		m = mdb.models['Model-1']
		disp_top = float(disp_top)
		m.boundaryConditions['BC-top'].setValues(u2 = disp_top)

		j = mdb.Job(activateLoadBalancing=False, atTime=None, contactPrint=OFF, 
			description='', echoPrint=OFF, explicitPrecision=DOUBLE_PLUS_PACK, 
			historyPrint=OFF, memory=90, memoryUnits=PERCENTAGE, model='Model-1', 
			modelPrint=OFF, multiprocessingMode=DEFAULT, name = jname, 
			nodalOutputPrecision=SINGLE, numCpus=num_threads, numDomains=num_threads, 
			parallelizationMethodExplicit=DOMAIN, queue=None, resultsFormat=ODB, 
			scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

		j.writeInput()

	#end
#end;