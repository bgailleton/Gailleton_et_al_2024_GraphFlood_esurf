'''
High level grpahflood object
'''
import dagger as dag
import scabbard as scb
import numpy as np
import matplotlib.pyplot as plt
import cmcrameri.cm as cm




class GraphFlood(object):

	"""
		docstring for GraphFlood
	"""

	def __init__(self, 
		grid,
		verbose = False,
		convergence_tracking = True,
		**kwargs
		):

		super(GraphFlood, self).__init__()
		
		# Grid object
		self.grid = grid

		# Checking and or feeding the graph
		if (self.grid.con is None and self.grid.graph is None):
			self.grid.compute_graphcon()
		elif (self.grid.con is not None and self.grid.graph is None):
			self.grid.graph = dag.graph(self.grid.con)
			self.grid.graph.compute_graph(self._Z, False, True)
		
		# Initialising the c++ graphflood object
		self.flood = dag.graphflood(self.grid.graph, self.grid.con)

		# feeding the topo
		self.flood.set_topo(self.grid._Z)

		# run all the configuration options
		self.config(**kwargs)

		self.active_figs = []
		self._callbacks = []

		self.hydro_dt = 1e-3 	

		self.verbose = verbose

		self.update_grid = True


		# time
		self.cumulative_time_hydro = 0.
		self.nit_hydro = 0
		self.cumulative_time_morpho = 0.

		# Convergence criteriae
		self.init_convergence_tracking(convergence_tracking)

	def config(self, **kwargs):
		'''
			Configures some of the general model options and inputs.
			Wrapper on the c++ object, not all options are in to keep things clear.
			Options:
				- dt or hydro_dt: hydrologic time ste
				- SFD (bool): if True set the mode to single flow, else multiple flow
				- minima (str): configure the method used for local minima resolution
					+ "reroute": reroute flow without filling the minimas
					+ "ignore": stop the flow at local minimas
					+ "fill" (default AND reccomended): fill the local minima with water

		'''

		if("SFD" in kwargs.keys()):
			if(kwargs["SFD"] == True):
				print("set flow conditions to single flow direction") if self.verbose else 0
				self.flood.enable_SFD()
			else:
				print("set flow conditions to mutliple flow direction") if self.verbose else 0
				self.flood.enable_MFD()

		if("minima" in kwargs.keys()):
			if(kwargs["minima"].lower() == "reroute"):
				self.flood.reroute_minima()
			elif(kwargs["minima"].lower() == "ignore"):
				self.flood.ignore_minima()
			else:
				self.flood.fill_minima()

		if("dt" in kwargs.keys()):
			self.flood.set_dt_hydro(kwargs["dt"])
			self.hydro_dt = kwargs["dt"]

		if("hydro_dt" in kwargs.keys()):
			self.flood.set_dt_hydro(kwargs["hydro_dt"])
			self.hydro_dt = kwargs["hydro_dt"]

		if("n_trackers" in kwargs.keys()):
			self.n_trackers = kwargs['n_trackers']

	def init_convergence_tracking(self, yes = True):

		self.convergence_trackers = yes
		self.n_trackers = 50
		self.n_pits = []
		self.dhw_monitoring = []
		self.Qwratio = []

		if(yes):
			self.flood.enable_Qwout_recording()
			self.flood.enable_dhw_recording()

	def set_precipitations(self,values = 1e-4):

		if(isinstance(values, np.ndarray)):
			if(np.prod(values.shape) == self.grid.nxy ):
				self.flood.set_water_input_by_variable_precipitation_rate(values.ravel())
			else:
				raise RuntimeError("array of precipitations needs to be of grid size or a single scalar")

		else:
			self.flood.set_water_input_by_constant_precipitation_rate(values)


	def set_input_discharge(self, nodes, values, asprec = False):
		'''
			Sets the input discharge at given locations (nodes need to be in flattened node indices)
		'''

		if(isinstance(nodes,list)):
			nodes = np.array(nodes, dtype=np.int32)
		if(isinstance(values,list)):
			values = np.array(values, dtype=np.float64)
		if(np.prod(nodes.shape) != np.prod(values.shape)):
			raise RuntimeError(f"nodes and values need to be the same size. Currently {np.prod(nodes)} vs {np.prod(values)}")

		# YOLO
		if(asprec == False):
			self.flood.set_water_input_by_entry_points(values, nodes)
		else:
			P = np.zeros(self.grid.nxy)
			P[nodes] = values
			self.flood.set_water_input_by_variable_precipitation_rate( P.ravel())


	def run_hydro(self, n_steps = 1, fig_update_step = 1,
	 force_morpho = False, run_morpho_every = 5, min_iteration_morpho = 500, run_morpho_ratio = 0.95,
	 runtime_callback = [], RAT = False, RAT_step = 1000, courant = True, **kwargs):
		'''
		'''
		# Only hydro on this version
		if(force_morpho):
			self.flood.enable_morpho()
		else:
			self.flood.disable_morpho()

		if(courant):
			self.flood.enable_courant_dt_hydro()
		else:
			self.flood.set_dt_hydro(self.hydro_dt)
		# print("DEBUG::DT is", self.hydro_dt)

		# Running loop
		for i in range(n_steps):

			# if i> min_iteration_morpho and i%run_morpho_every == 0 and force_morpho:
			totin = self.flood.get_tot_Qw_input()
			if(totin >0 and self.flood.get_tot_Qw_output() / totin > run_morpho_ratio and force_morpho):
				# print("DEGUB::MORPHO ON")
				self.flood.enable_morpho()
				if("block_uplift_rate" in kwargs):
					self.flood.block_uplift(kwargs["block_uplift_rate"])
			else:
				self.flood.disable_morpho()
			# if(i > 5000 and i%5 == 0):
			# 	input("tttt")

			self.cumulative_time_hydro += self.flood.get_dt_hydro()
			self.nit_hydro += 1

			# Running hte actual model
			self.flood.run()

			# Balance checker for debugging purposes
			balance = abs(self.flood.get_tot_Qw_input() - self.flood.get_tot_Qwin_output())
			if(balance > 1):
				print("WARNING::Qw-in imbalance -> " + str(self.flood.get_tot_Qw_input()) + " got in and " + str(self.flood.get_tot_Qwin_output()) + " got out. Unbalance is: " + str(balance))

			if(i%500 == 0):
				print(i)

			if(RAT):
				if(i>0 and i%RAT_step == 0):
					newdt = input("new_dt: ")
					if(newdt == ''):
						newdt = self.hydro_dt
					else:
						self.hydro_dt = float(newdt)

					self.flood.set_dt_hydro(self.hydro_dt)


			for rtcb in runtime_callback:
				if(isinstance(rtcb,list)):
					rtcb[0](*rtcb[1])
				else:
					rtcb()

			# Monitoring the water convergence now:
			if(self.convergence_trackers):
				if(len(self.n_pits) == self.n_trackers):
					self.n_pits.pop(0)
					self.dhw_monitoring.pop(0)
					self.Qwratio.pop(0)

				self.n_pits.append(self.grid.graph.get_n_pits())
				arr = self.flood.get_dhw_recording()
				self.dhw_monitoring.append([np.percentile(arr,5), np.percentile(arr,25), np.median(arr), np.percentile(arr,75), np.percentile(arr,95)])
				arr = self.flood.get_Qwout_recording()
				mask = arr == 0
				arr = self.flood.get_Qwin()/arr
				arr[mask] = 1
				self.Qwratio.append([np.percentile(arr,5), np.percentile(arr,25), np.median(arr), np.percentile(arr,75), np.percentile(arr,95)])

				if(self.verbose):
					print("\n###############################")
					print("###############################")
					print("Monitoring results:")
					print("N pits:",self.n_pits[-1])
					print("delta hw:",self.dhw_monitoring[-1])
					print("Qwratio:",self.Qwratio[-1])
					print("###############################")
					print("###############################\n")

			

			# updating the figures
			if(i%fig_update_step == 0):
				self.update_figs()
				if self.update_grid:
					self.grid._Z = self.flood.get_bedrock_topo()
		if self.update_grid:
			self.grid._Z = self.flood.get_bedrock_topo()

	def update_figs(self):
		for tax in self._callbacks:
			tax.update()
		for tf in self.active_figs:
			tf.canvas.draw_idle()
			tf.canvas.start_event_loop(0.001)


	def get_xmonitor_HYDRO(self):
		tx1 = max(self.nit_hydro - self.n_trackers, 0)
		tx2 = self.nit_hydro
		# print("DEBUG_XMONO::",tx1,tx2)
		return np.arange(tx1 + 1, tx2 + 1)

	def get_monitor_pits(self):
		return [self.get_xmonitor_HYDRO(), self.n_pits]

	def get_monitor_median_hw_pits(self):
		# print(self.dhw_monitoring[2], "KKKK")
		return [self.get_xmonitor_HYDRO(), np.array(self.dhw_monitoring)[:,2]]

	def get_monitor_10th_hw_pits(self):
		return [self.get_xmonitor_HYDRO(), np.array(self.dhw_monitoring)[:,0]]


	def get_monitor_90th_hw_pits(self):
		return [self.get_xmonitor_HYDRO(), np.array(self.dhw_monitoring)[:,4]]

	def debugyolo(self):
		return np.array(self.grid.graph.get_debug_mask())

	def pop_Qw_fig(self, jupyter = False, clim = None, alpha_hillshade = 0.3):

		if(jupyter == False):
			plt.ioff()

		fig,ax = plt.subplots()
		dax = scb.RGridDax(self.grid, ax, alpha_hillshade=1)
		arr = self.flood.get_Qwin()
		if(np.prod(arr.shape) == 0):
			arr = self.grid.zeros()
		else:
			arr.reshape(self.grid.rshp)

		Qwax = dax.drape_on(arr, cmap = "Blues", clim = clim, delta_zorder = 1, alpha = 1 - alpha_hillshade, callback = self.flood.get_Qwin)
		# Qwax = dax.drape_on(arr, cmap = "Blues", clim = None, delta_zorder = 1, alpha = 0.9, callback = self.debugyolo)

		self.active_figs.append(fig)
		self._callbacks.append(dax)
		self._callbacks.append(Qwax)
		fig.show()
		fig.canvas.draw_idle()
		fig.canvas.start_event_loop(0.001)

	def pop_hw_fig(self, jupyter = False, clim = None, alpha_hillshade = 0.3):

		if(jupyter == False):
			plt.ioff()

		fig,ax = plt.subplots()
		dax = scb.RGridDax(self.grid, ax, alpha_hillshade=1)
		arr = self.flood.get_Qwin()
		if(np.prod(arr.shape) == 0):
			arr = self.grid.zeros()
		else:
			arr.reshape(self.grid.rshp)

		Qwax = dax.drape_on(arr, cmap = "Blues", clim = clim, delta_zorder = 1, alpha = 1 - alpha_hillshade, callback = self.flood.get_hw)

		plt.colorbar(Qwax.im, label = "Water depth (m)")
		self.active_figs.append(fig)
		self._callbacks.append(dax)
		self._callbacks.append(Qwax)
		fig.show()
		fig.canvas.draw_idle()
		fig.canvas.start_event_loop(0.001)


	def pop_custom_fig(self, jupyter = False, clim = None, callback = None,callback_params = None, cmap = 'magma', alpha_hillshade = 0.3):

		if(callback is None):
			raise RuntimeError("Need callback for custom fig")

		if(jupyter == False):
			plt.ioff()

		fig,ax = plt.subplots()
		dax = scb.RGridDax(self.grid, ax, alpha_hillshade=1)
		arr = self.flood.get_Qwin()
		if(np.prod(arr.shape) == 0):
			arr = self.grid.zeros()
		else:
			arr.reshape(self.grid.rshp)

		Qwax = dax.drape_on(arr, cmap = cmap, clim = clim, delta_zorder = 1, alpha = 1 - alpha_hillshade, callback = callback, callback_params = callback_params)

		plt.colorbar(Qwax.im, label = "Water depth (m)")
		self.active_figs.append(fig)
		self._callbacks.append(dax)
		self._callbacks.append(Qwax)
		fig.show()
		fig.canvas.draw_idle()
		fig.canvas.start_event_loop(0.001)


	def pop_topo_fig(self, jupyter = False, clim = None, bedrock = False, alpha_hillshade =0.5):

		if(jupyter == False):
			plt.ioff()

		fig,ax = plt.subplots()
		dax = scb.RGridDax(self.grid, ax, alpha_hillshade=1)
		arr = self.flood.get_Qwin()
		if(np.prod(arr.shape) == 0):
			arr = self.grid.zeros()
		else:
			arr.reshape(self.grid.rshp)

		topo = dax.drape_on(arr, cmap = "gist_earth", clim = clim, delta_zorder = 1, alpha =1 - alpha_hillshade, callback = self.flood.get_bedrock_topo if bedrock else self.flood.get_surface_topo)

		plt.colorbar(topo.im, label = "Surface (hw + Z) (m)" if bedrock == False else "Z (m)")
		self.active_figs.append(fig)
		self._callbacks.append(dax)
		self._callbacks.append(topo)
		fig.show()
		fig.canvas.draw_idle()
		fig.canvas.start_event_loop(0.001)

	def pop_dtopo_fig(self, jupyter = False, clim = None, alpha_hillshade = 0.5):

		if(jupyter == False):
			plt.ioff()

		self.otopo = np.copy(self.grid.Z2D)

		fig,ax = plt.subplots()
		dax = scb.RGridDax(self.grid, ax, alpha_hillshade=1)
		arr = self.flood.get_Qwin()
		if(np.prod(arr.shape) == 0):
			arr = self.grid.zeros()
		else:
			arr.reshape(self.grid.rshp)

		def _pop_dtopo_fig():
			return self.otopo - self.flood.get_bedrock_topo().reshape(self.grid.rshp)

		topo = dax.drape_on(arr, cmap = cm.bam, clim = clim, delta_zorder = 1, alpha = 1 - alpha_hillshade, callback = _pop_dtopo_fig)

		plt.colorbar(topo.im, label = "$\Delta Z$")
		self.active_figs.append(fig)
		self._callbacks.append(dax)
		self._callbacks.append(topo)
		fig.show()
		fig.canvas.draw_idle()
		fig.canvas.start_event_loop(0.001)



	def pop_monitor_fig(self, jupyter = False):

		if(self.convergence_trackers == False):
			raise AttributeError("cannot pop monitor figure if convergence_trackers is not activated")

		if(jupyter == False):
			plt.ioff()

		# # gridspec inside gridspec
		# fig = plt.figure()
		# gs0 = gridspec.GridSpec(1, 2, figure=fig)

		fig,ax1 = plt.subplots()

		ax2 = ax1.twinx()

		ax1.set_xlabel('N iterations')
		ax1.set_ylabel('N internal pits')

		l1 = ax1.plot([0],[0], lw = 2, color = 'r')

		l3 = ax2.plot([0],[0], lw = 1, ls = '--', color = 'k')
		l4 = ax2.plot([0],[0], lw = 1, ls = '--', color = 'k')
		l2 = ax2.plot([0],[0], lw = 2, color = 'k')

		dax1 = scb.callbax_sline(ax1, l1, self.get_monitor_pits, axylim = None)

		dax2 = scb.callbax_sline(ax2, l2, self.get_monitor_median_hw_pits, axylim = None)
		dax3 = scb.callbax_sline(ax2, l3, self.get_monitor_10th_hw_pits, axylim = None, axylim_ignore = True)
		dax4 = scb.callbax_sline(ax2, l4, self.get_monitor_90th_hw_pits, axylim = None, axylim_ignore = True)
		


		self.active_figs.append(fig)
		self._callbacks.append(dax1)
		self._callbacks.append(dax2)
		self._callbacks.append(dax3)
		self._callbacks.append(dax4)
		# self._callbacks.append(topo)
		fig.show()
		fig.canvas.draw_idle()
		fig.canvas.start_event_loop(0.001)


	def pop_Qwratio_fig(self, jupyter = False, clim = None, alpha_hillshade = 0.3):

		self.flood.enable_Qwout_recording()

		if(jupyter == False):
			plt.ioff()

		fig,ax = plt.subplots()
		dax = scb.RGridDax(self.grid, ax, alpha_hillshade=1)
		arr = self.flood.get_Qwin()
		if(np.prod(arr.shape) == 0):
			arr = self.grid.zeros()
		else:
			arr.reshape(self.grid.rshp)

		def _cback():
			Qwin = self.flood.get_Qwin()
			mask = Qwin == 0
			tarr = self.flood.get_Qwout_recording()/Qwin
			tarr[mask] = 1.
			return tarr

		Qwax = dax.drape_on(arr, cmap = "RdBu_r", clim = (0.75,1.25), delta_zorder = 1, alpha = 1 - alpha_hillshade, callback = _cback)
		# Qwax = dax.drape_on(arr, cmap = "Blues", clim = None, delta_zorder = 1, alpha = 0.9, callback = self.debugyolo)

		plt.colorbar(Qwax.im, label = "$Q_{wout}/Q_{win}")

		self.active_figs.append(fig)
		self._callbacks.append(dax)
		self._callbacks.append(Qwax)
		fig.show()
		fig.canvas.draw_idle()
		fig.canvas.start_event_loop(0.001)






















































#End of file