"""
# 创建约束，指定Q值范围
q_range = {
    'qmin': 0.0001,  # Å^-1
    'qmax': 30.0,    # Å^-1
    'qstep': 0.01    # Å^-1
}

# 创建PDF约束并添加到引擎
pdc = PairDistributionConstraint(
    experimentalData="g_r.dat",  # 使用g(r)数据而非G(r)数据
    weighting="xrays_scattering",  # 使用X射线散射因子
    qRange=q_range      # 指定Q值范围

修改取消Gr，使用gr作为约束，函数Gr不变

PairDistributionConstraints contains classes for all constraints related
to experimental pair distribution functions.

.. inheritance-diagram:: fullrmc.Constraints.PairDistributionConstraints
    :parts: 1
"""
# standard libraries imports
from __future__ import print_function
import itertools, inspect, copy, os, re

# external libraries imports
import numpy as np
from pdbparser.Utilities.Database import is_element_property, get_element_property
from pdbparser.Utilities.Collection import get_normalized_weighting

# fullrmc imports
from ..Globals import INT_TYPE, FLOAT_TYPE, PI, PRECISION, LOGGER
from ..Globals import str, long, unicode, bytes, basestring, range, xrange, maxint
from ..Core.Collection import is_number, is_integer, get_path
from ..Core.Collection import reset_if_collected_out_of_date, get_real_elements_weight
from ..Core.Collection import get_caller_frames
from ..Core.Constraint import Constraint, ExperimentalConstraint
from ..Core.pairs_histograms import multiple_pairs_histograms_coords, full_pairs_histograms_coords
from ..Constraints.Collection import ShapeFunction


class PairDistributionConstraint(ExperimentalConstraint):
    """
    控制原子配置的对分布函数(PDF)，表示为g(r)。对分布函数（又称径向分布函数）
    直接从粉末衍射实验数据中计算获得。
    
    g(r)是通过将所有原子对距离进行加权直方图统计来计算的，得到局部数密度，
    具体计算如下：

    .. math::
        g(r) = \\sum \\limits_{i,j}^{N} w_{i,j} \\frac{\\rho_{i,j}(r)}{\\rho_{0}}
             = \\sum \\limits_{i,j}^{N} w_{i,j} \\frac{n_{i,j}(r) / v(r)}{N_{i,j} / V}

    其中：
    :math:`r` 是两原子之间的距离。
    :math:`\\rho_{i,j}(r)` 是原子i和j的对密度函数。
    :math:`\\rho_{0}` 是系统的平均数密度。
    :math:`w_{i,j}` 是原子类型i和j的相对权重。
    :math:`N` 是原子总数。
    :math:`V` 是系统体积。
    :math:`n_{i,j}(r)` 是与j原子距离为r的i原子数量。
    :math:`v(r)` 是距离r处厚度为dr的球壳体积。
    :math:`N_{i,j}` 是系统中原子i和j的总数量。
    
    注意：此版本的约束直接使用g(r)作为约束数据，而不是转换为G(r)。因此，
    实验输入数据也应该是g(r)格式，而不是G(r)格式。

    :Parameters:
        #. experimentalData (numpy.ndarray, string): Experimental g(r) data as
           numpy.ndarray or string path to load data using numpy.loadtxt.
        #. dataWeights (None, numpy.ndarray): Weights array of the same number
           of points of experimentalData used in the constraint's standard
           error computation. Therefore particular fitting emphasis can be
           put on different data points that might be considered as more or less
           important in order to get a reasonable and plausible modal.\n
           If None is given, all data points are considered of the same
           importance in the computation of the constraint's standard error.\n
           If numpy.ndarray is given, all weights must be positive and all
           zeros weighted data points won't contribute to the total
           constraint's standard error. At least a single weight point is
           required to be non-zeros and the weights array will be automatically
           scaled upon setting such as the the sum of all the weights
           is equal to the number of data points.
        #. weighting (string): The elements weighting scheme. It must be any
           atomic attribute (atomicNumber, neutronCohb, neutronIncohb,
           neutronCohXs, neutronIncohXs, atomicWeight, covalentRadius) defined
           in pdbparser database. In case of xrays or neutrons experimental
           weights, one can simply set weighting to 'xrays' or 'neutrons'
           and the value will be automatically adjusted to respectively
           'atomicNumber' and 'neutronCohb'. If attribute values are
           missing in the pdbparser database, atomic weights must be
           given in atomsWeight dictionary argument.
        #. atomsWeight (None, dict): Atoms weight dictionary where keys are
           atoms element and values are custom weights. If None is given
           or partially given, missing elements weighting will be fully set
           given weighting scheme.
        #. scaleFactor (number): A normalization scale factor used to normalize
           the computed data to the experimental ones.
        #. adjustScaleFactor (list, tuple): Used to adjust fit or guess
           the best scale factor during stochastic engine runtime.
           It must be a list of exactly three entries.\n
           #. The frequency in number of generated moves of finding the best
              scale factor. If 0 frequency is given, it means that the scale
              factor is fixed.
           #. The minimum allowed scale factor value.
           #. The maximum allowed scale factor value.
        #. shapeFuncParams (None, numpy.ndarray, dict): The shape function is
           subtracted from the total g(r). It must be used when non-periodic
           boundary conditions are used to take into account the atomic
           density drop and to correct for the :math:`\\rho_{0}` approximation.
           The shape function can be set to None which means unsused, or set
           as a constant shape given by a numpy.ndarray or computed
           from all atoms and updated every 'updateFreq' accepted moves.
           If dict is given the following keywords can be given, otherwise
           default values will be automatically set.\n
           * **rmin (number) default (0.00) :** The minimum distance in
             :math:`\\AA` considered upon building the histogram prior to
             computing the shape function. If not defined, rmin will be
             set automatically to 0.
           * **rmax (None, number) default (None) :** The maximum distance
             in :math:`\\AA` considered upon building the histogram prior
             to computing the shape function. If not defnined, rmax will
             be automatically set to :math:`maximum\ box\ length + 10\\AA`
             at engine runtime.
           * **dr (number) default (0.5) :** The bin size in :math:`\\AA`
             considered upon building the histogram prior to computing the
             shape function. If not defined, it will be automatically set
             to 0.5.
           * **qmin (number) default (0.001) :** The minimum reciprocal
             distance q in :math:`\\AA^{-1}` considered to compute the
             shape function. If not defined, it will be automatically
             set to 0.001.
           * **qmax (number) default (0.75) :** The maximum reciprocal
             distance q in :math:`\\AA^{-1}` considered to compute the
             shape function. If not defined, it will be automatically
             set to 0.75.
           * **dq (number) default (0.005) :** The reciprocal distance bin
             size in :math:`\\AA^{-1}` considered to compute the shape
             function. If not defined, it will be automatically
             set to 0.005.
           * **updateFreq (integer) default (1000) :** The frequency of
             recomputing the shape function in number of accpeted moves.
        #. windowFunction (None, numpy.ndarray): The window function to
           convolute with the computed pair distribution function of the
           system prior to comparing it with the experimental data.
        #. limits (None, tuple, list): The distance limits to compute the
           histograms. If None is given, the limits will be automatically
           set the the min and max distance of the experimental data.
           Otherwise, a tuple of exactly two items where the first is the
           minimum distance or None and the second is the maximum distance
           or None.

    **NB**: If adjustScaleFactor first item (frequency) is 0, the scale factor
    will remain untouched and the limits minimum and maximum won't be checked.

    .. code-block:: python

        # import fullrmc modules
        from fullrmc.Engine import Engine
        from fullrmc.Constraints.PairDistributionConstraints import PairDistributionConstraint

        # create engine
        ENGINE = Engine(path='my_engine.rmc')

        # set pdb file
        ENGINE.set_pdb('system.pdb')

        # create and add constraint
        PDC = PairDistributionConstraint(experimentalData="g_r.dat", weighting="atomicNumber")
        ENGINE.add_constraints(PDC)
    """
    # Waasmaier-Kirfel参数 - 用于X射线散射因子计算
    _wk_params = {
        "C": ([2.657506, 1.078079, 1.490909, -4.241070, 0.713791],
              [14.780758, 0.776775, 42.086842, -0.000294, 0.239535],
              4.297983),
        "H": ([0.413048, 0.294953, 0.187491, 0.080701, 0.023736],
              [15.569946, 32.398468, 5.711404, 61.889874, 1.334118],
              0.000049),
        "O": ([2.960427, 2.508818, 0.637853, 0.722838, 1.142756],
              [14.182259, 5.936858, 0.112726, 34.958481, 0.390240],
              0.027014),
        "N": ([11.893780, 3.277479, 1.858092, 0.858927, 0.912985],
              [0.000158, 10.232723, 30.344690, 0.656065, 0.217287],
              -11.804902),
        "S": ([6.372157, 5.154568, 1.473732, 1.635073, 1.209372],
              [1.514347, 22.092527, 0.061373, 55.445175, 0.646925],
              0.154722)
        # 可以根据需要添加更多元素的参数
    }
    
    # Q值范围默认设置
    _q_range_defaults = {
        'qmin': 0.0001,  # Å^-1
        'qmax': 30.0,    # Å^-1
        'qstep': 0.01    # Å^-1
    }
    
    def __init__(self, experimentalData, dataWeights=None, weighting="atomicNumber",
                       atomsWeight=None, scaleFactor=1.0, adjustScaleFactor=(0, 0.8, 1.2),
                       shapeFuncParams=None, windowFunction=None, limits=None,qRange=None):
        # initialize constraint
        super(PairDistributionConstraint, self).__init__(experimentalData=experimentalData, dataWeights=dataWeights, scaleFactor=scaleFactor, adjustScaleFactor=adjustScaleFactor)

        # 设置Q值范围
        self._q_range = self._q_range_defaults.copy()
        if qRange is not None:
            self.set_q_range(qRange)

        # set elements weighting
        self.set_weighting(weighting)
        # set atomsWeight
        self.set_atoms_weight(atomsWeight)
        # set window function
        self.set_window_function(windowFunction)
        # set shape function parameters
        self.set_shape_function_parameters(shapeFuncParams)

        # set frame data
        FRAME_DATA = [d for d in self.FRAME_DATA]
        FRAME_DATA.extend(['_PairDistributionConstraint__bin',
                           '_PairDistributionConstraint__minimumDistance',
                           '_PairDistributionConstraint__maximumDistance',
                           '_PairDistributionConstraint__shellCenters',
                           '_PairDistributionConstraint__edges',
                           '_PairDistributionConstraint__histogramSize',
                           '_PairDistributionConstraint__experimentalDistances',
                           '_PairDistributionConstraint__experimentalPDF',
                           '_PairDistributionConstraint__shellVolumes',
                           '_PairDistributionConstraint__elementsPairs',
                           '_PairDistributionConstraint__weighting',
                           '_PairDistributionConstraint__atomsWeight',
                           '_PairDistributionConstraint__weightingScheme',
                           '_PairDistributionConstraint__windowFunction',
                           '_q_range',
                           '_elementsWeight',
                           '_shapeFuncParams',
                           '_shapeUpdateFreq',
                           '_shapeArray', ] )

        RUNTIME_DATA = [d for d in self.RUNTIME_DATA]
        RUNTIME_DATA.extend( ['_shapeArray'] )
        object.__setattr__(self, 'FRAME_DATA',   tuple(FRAME_DATA)  )
        object.__setattr__(self, 'RUNTIME_DATA', tuple(RUNTIME_DATA) )

    def set_q_range(self, qRange):
        """
        设置用于散射因子计算的Q值范围。
    
        :Parameters:
            #. qRange (dict): 包含Q值范围设置的字典，可以包含以下键：
               * qmin (float): 最小Q值，单位为Å^-1
               * qmax (float): 最大Q值，单位为Å^-1
               * qstep (float): Q值步长，单位为Å^-1
        """
        assert isinstance(qRange, dict), LOGGER.error("qRange must be a dictionary")
    
        for key in qRange:
            assert key in self._q_range, LOGGER.error(f"Unknown Q range parameter: {key}")
            try:
                value = float(qRange[key])
                assert value > 0, LOGGER.error(f"Q range parameter {key} must be positive")
                self._q_range[key] = value
            except:
                raise LOGGER.error(f"Q range parameter {key} must be a positive number")
    
        # 验证qmin < qmax
        assert self._q_range['qmin'] < self._q_range['qmax'], LOGGER.error("qmin must be less than qmax")
    
        # 验证步长合理
        num_steps = (self._q_range['qmax'] - self._q_range['qmin']) / self._q_range['qstep']
        assert num_steps <= 10000, LOGGER.error("Too many Q points. Consider increasing qstep or reducing qmax")
    
        # 转储到存储库
        self._dump_to_repository({'_PairDistributionConstraint__q_range': self._q_range})

    def _calculate_scattering_factors(self, elements=None):
        """
        计算给定元素在Q值范围内的散射因子。
    
        :Parameters:
            #. elements (list, None): 要计算散射因子的元素列表。如果为None，则使用engine.elements。
    
        :Returns:
            #. scatter_factors (dict): 散射因子字典，键为元素，值为散射因子数组
        """
        if elements is None:
            elements = self.engine.elements
    
        # 准备Q值数组
        qmin = self._q_range['qmin']
        qmax = self._q_range['qmax']
        qstep = self._q_range['qstep']
        num_q_points = int((qmax - qmin) / qstep) + 1
        q_values = np.linspace(qmin, qmax, num_q_points)
    
        # 初始化散射因子字典
        scatter_factors = {}
    
        # 计算每个元素在每个Q值处的散射因子
        for element in elements:
            # 初始化散射因子数组
            f_q = np.zeros(num_q_points, dtype=FLOAT_TYPE)
        
            if element in self._wk_params:
                # 获取Waasmaier-Kirfel参数
                a_vals, b_vals, c_val = self._wk_params[element]
            
                # 对每个Q值计算散射因子
                for i, q in enumerate(q_values):
                    # 计算s^2，其中s = q/(4π)
                    s2 = (q * q) / (16.0 * PI * PI)
                
                    # 使用Waasmaier-Kirfel公式计算散射因子
                    f_q[i] = c_val
                    for k in range(5):
                        f_q[i] += a_vals[k] * np.exp(-b_vals[k] * s2)
            else:
                # 如果元素没有参数，使用原子序数作为近似
                try:
                    atomic_number = get_element_property(element=element, property="atomicNumber")
                    # 对于缺少参数的元素，假设散射因子简单地随Q衰减
                    f_q = np.ones(num_q_points, dtype=FLOAT_TYPE) * atomic_number
                    # 简单衰减模型
                    f_q *= np.exp(-q_values/20.0)
                except:
                    LOGGER.warn(f"元素 {element} 未找到参数，使用1.0作为散射因子")
                    f_q = np.ones(num_q_points, dtype=FLOAT_TYPE)
        
            # 存储散射因子
            scatter_factors[element] = f_q
    
        return scatter_factors, q_values

    def _calculate_xray_weights(self):
        """计算X射线散射权重
    
        使用Waasmaier-Kirfel参数化计算原子散射因子
        先计算每个Q值点的平均散射因子
        然后在Q空间进行积分平均得到最终权重
        """
        elements = self.engine.elements
        fractions = {}  # 原子分数字典
    
        # 计算原子分数
        total_atoms = sum(self.engine.numberOfAtomsPerElement.values())
        for elem in elements:
            fractions[elem] = self.engine.numberOfAtomsPerElement[elem] / total_atoms
    
        # 计算散射因子
        scatter_factors, q_values = self._calculate_scattering_factors(elements)
        num_q_points = len(q_values)
       
        # 第一步：计算每个Q值点的平均散射因子
        f_av_q = np.zeros(num_q_points, dtype=FLOAT_TYPE)
        for q_idx in range(num_q_points):
            for elem in elements:
                f_i = scatter_factors[elem][q_idx] 
                f_av_q[q_idx] += fractions[elem] * f_i
    
        # 计算平均散射因子的平方（用于归一化）
        f_av_q_sq = f_av_q * f_av_q
    
        # 初始化原子对权重字典（存储每个Q值的权重）
        # 键格式：'elem_i-elem_j'，值：每个Q点的权重数组
        pair_weights_q = {}
        for i, elem_i in enumerate(elements):
            for j in range(i, len(elements)):
                elem_j = elements[j]
                pair_key = f"{elem_i}-{elem_j}"
                pair_weights_q[pair_key] = np.zeros(num_q_points, dtype=FLOAT_TYPE)
    
        # 第二步：计算每对原子在每个Q值点的权重
        for q_idx in range(num_q_points):
            # 遍历所有原子对
            for i, elem_i in enumerate(elements):
                f_i = scatter_factors[elem_i][q_idx] 
            
                for j in range(i, len(elements)):
                    elem_j = elements[j]
                    f_j = scatter_factors[elem_j][q_idx] 
                
                    # 计算原子对权重
                    pair_key = f"{elem_i}-{elem_j}"
                    weight = f_i * f_j * fractions[elem_i] * fractions[elem_j]
                
                    # 归一化
                    if f_av_q_sq[q_idx] > 0:  # 避免除零错误
                        weight /= f_av_q_sq[q_idx]
                
                    # 对不同类型原子对，权重乘以2
                    if elem_i != elem_j:
                        weight *= 2.0
                
                    # 存储该Q值点的权重
                    pair_weights_q[pair_key][q_idx] = weight
    
        # 第三步：对Q空间进行积分平均得到最终权重
        weightingScheme = {}
    
        for pair_key, weights in pair_weights_q.items():
            # 使用梯形积分法在Q空间积分
            total_weight = 0.0
            for q_idx in range(num_q_points - 1):
                dq = q_values[q_idx + 1] - q_values[q_idx]
                w1 = weights[q_idx]
                w2 = weights[q_idx + 1]
                total_weight += 0.5 * dq * (w1 + w2)
        
            # 归一化积分结果
            q_range = q_values[-1] - q_values[0]
            if q_range > 0:  # 避免除零错误
                total_weight /= q_range
        
            # 存储最终权重
            weightingScheme[pair_key] = FLOAT_TYPE(total_weight)
    
        # 返回计算出的权重方案
        return weightingScheme

    def _codify_update__(self, name='constraint', addDependencies=True):
        dependencies = []
        code         = []
        if addDependencies:
            code.extend(dependencies)
        dw = self.dataWeights
        if dw is not None:
            dw = list(dw)
        code.append("dw = {dw}".format(dw=dw))
        sfp = self._shapeFuncParams
        if isinstance(sfp, np.ndarray):
            sfp = list(sfp)
        code.append("sfp = {sfp}".format(sfp=sfp))
        wf = self.windowFunction
        if isinstance(wf, np.ndarray):
            code.append("wf = np.array({wf})".format(wf=list(wf)))
        else:
            code.append("wf = {wf}".format(wf=wf))
        code.append("{name}.set_used({val})".format(name=name, val=self.used))
        code.append("{name}.set_scale_factor({val})".format(name=name, val=self.scaleFactor))
        code.append("{name}.set_adjust_scale_factor({val})".format(name=name, val=self.adjustScaleFactor))
        code.append("{name}.set_data_weights(dw)".format(name=name))
        code.append("{name}.set_atoms_weight({val})".format(name=name, val=self.atomsWeight))
        code.append("{name}.set_window_function(wf)".format(name=name))
        code.append("{name}.set_shape_function_parameters(sfp)".format(name=name))
        code.append("{name}.set_limits({val})".format(name=name, val=self.limits))
        # return
        return dependencies, '\n'.join(code)


    def _codify__(self, engine, name='constraint', addDependencies=True):
        assert isinstance(name, basestring), LOGGER.error("name must be a string")
        assert re.match('[a-zA-Z_][a-zA-Z0-9_]*$', name) is not None, LOGGER.error("given name '%s' can't be used as a variable name"%name)
        klass        = self.__class__.__name__
        dependencies = ['import numpy as np','from fullrmc.Constraints import {klass}s'.format(klass=klass)]
        code         = []
        if addDependencies:
            code.extend(dependencies)
        x = list(self.experimentalData[:,0])
        y = list(self.experimentalData[:,1])
        code.append("x = {x}".format(x=x))
        code.append("y = {y}".format(y=y))
        code.append("d = np.transpose([x,y]).astype(np.float32)")
        dw = self.dataWeights
        if dw is not None:
            dw = list(dw)
        code.append("dw = {dw}".format(dw=dw))
        sfp = self._shapeFuncParams
        if isinstance(sfp, np.ndarray):
            sfp = list(sfp)
        code.append("sfp = {sfp}".format(sfp=sfp))
        wf = self.windowFunction
        if isinstance(wf, np.ndarray):
            code.append("wf = np.array({wf})".format(wf=list(wf)))
        else:
            code.append("wf = {wf}".format(wf=wf))
        code.append("{name} = {klass}s.{klass}\
    (experimentalData=d, dataWeights=dw, weighting='{weighting}', atomsWeight={atomsWeight}, \
    scaleFactor={scaleFactor}, adjustScaleFactor={adjustScaleFactor}, \
    shapeFuncParams=sfp, windowFunction=wf, limits={limits}, qRange={qRange})".format(name=name, klass=klass,
                    weighting=self.weighting, atomsWeight=self.atomsWeight,
                    scaleFactor=self.scaleFactor, adjustScaleFactor=self.adjustScaleFactor,
                    shapeFuncParams=sfp, limits=self.limits, qRange=self._q_range))
        code.append("{engine}.add_constraints([{name}])".format(engine=engine, name=name))
        # return
        return dependencies, '\n'.join(code)

    def _on_collector_reset(self):
        pass

    def _update_shape_array(self):
        rmin = self._shapeFuncParams['rmin']
        rmax = self._shapeFuncParams['rmax']
        dr   = self._shapeFuncParams['dr'  ]
        qmin = self._shapeFuncParams['qmin']
        qmax = self._shapeFuncParams['qmax']
        dq   = self._shapeFuncParams['dq'  ]
        if rmax is None:
            if self.engine.isPBC:
                a = self.engine.boundaryConditions.get_a()
                b = self.engine.boundaryConditions.get_b()
                c = self.engine.boundaryConditions.get_c()
                rmax = FLOAT_TYPE( np.max([a,b,c]) + 10 )
            else:
                coordsCenter = np.sum(self.engine.realCoordinates, axis=0)/self.engine.realCoordinates.shape[0]
                coordinates  = self.engine.realCoordinates-coordsCenter
                distances    = np.sqrt( np.sum(coordinates**2, axis=1) )
                maxDistance  = 2.*np.max(distances)
                rmax = FLOAT_TYPE( maxDistance + 10 )
                LOGGER.warn("@%s Better set shape function rmax with infinite boundary conditions. Here value is automatically set to %s"%(self.engine.usedFrame, rmax))
        shapeFunc  = ShapeFunction(engine    = self.engine,
                                   weighting = self.__weighting,
                                   qmin=qmin, qmax=qmax, dq=dq,
                                   rmin=rmin, rmax=rmax, dr=dr)
        self._shapeArray = shapeFunc.get_Gr_shape_function( self.shellCenters )
        # 修改形状函数，从G(r)转换为g(r)形式，以便直接用于g(r)约束
        rho0 = self.engine.numberDensity
        self._shapeArray = self._shapeArray / (4.0 * PI * self.__shellCenters * rho0) + 1
        
        del shapeFunc
        # dump to repository
        self._dump_to_repository({'_shapeArray': self._shapeArray})

    def _reset_standard_error(self):
        # recompute squared deviation
        if self.data is not None:
            totalPDF = self.__get_total_Gr(self.data, rho0=self.engine.numberDensity)
            self.set_standard_error(self.compute_standard_error(modelData = totalPDF))

    def _runtime_initialize(self):
        if self._shapeFuncParams is None:
            self._shapeArray = None
        elif isinstance(self._shapeFuncParams, np.ndarray):
            self._shapeArray = self._shapeFuncParams[self.limitsIndexStart:self.limitsIndexEnd+1]
        elif isinstance(self._shapeFuncParams, dict) and self._shapeArray is None:
            self._update_shape_array()
        # reset standard error
        self._reset_standard_error()
        # set last shape update flag
        self._lastShapeUpdate = self.engine.accepted

    def _runtime_on_step(self):
        """ Update shape function when needed. and update engine total """
        if self._shapeUpdateFreq and self._shapeFuncParams is not None:
            if (self._lastShapeUpdate != self.engine.accepted) and not (self.engine.accepted%self._shapeUpdateFreq):
                # reset shape array
                self._update_shape_array()
                # reset standard error
                self._reset_standard_error()
                # update engine chiSquare
                oldTotalStandardError = self.engine.totalStandardError
                self.engine.update_total_standard_error()
                LOGGER.info("@%s Constraint '%s' shape function updated, engine total standard error updated from %.6f to %.6f" %(self.engine.usedFrame, self.__class__.__name__, oldTotalStandardError, self.engine.totalStandardError))
                self._lastShapeUpdate = self.engine.accepted

    @property
    def bin(self):
        """ Experimental data distances bin."""
        return self.__bin

    @property
    def minimumDistance(self):
        """ Experimental data minimum distances."""
        return self.__minimumDistance

    @property
    def maximumDistance(self):
        """ Experimental data maximum distances."""
        return self.__maximumDistance

    @property
    def histogramSize(self):
        """ Histogram size."""
        return self.__histogramSize

    @property
    def experimentalDistances(self):
        """ Experimental distances array."""
        return self.__experimentalDistances

    @property
    def shellCenters(self):
        """ Shells center array."""
        return self.__shellCenters

    @property
    def shellVolumes(self):
        """ Shells volume array."""
        return self.__shellVolumes

    @property
    def experimentalPDF(self):
        """ Experimental pair distribution function data."""
        return self.__experimentalPDF

    @property
    def elementsPairs(self):
        """ Elements pairs."""
        return self.__elementsPairs

    @property
    def weighting(self):
        """ Elements weighting definition."""
        return self.__weighting

    @property
    def atomsWeight(self):
        """Custom atoms weight"""
        return self.__atomsWeight

    @property
    def weightingScheme(self):
        """ Elements weighting scheme."""
        return self.__weightingScheme

    @property
    def windowFunction(self):
        """ Window function."""
        return self.__windowFunction

    @property
    def shapeArray(self):
        """ Shape function data array."""
        return self._shapeArray

    @property
    def shapeUpdateFreq(self):
        """Shape function update frequency."""
        return self._shapeUpdateFreq

    def _set_weighting_scheme(self, weightingScheme):
        """To be only used internally by PairDistributionConstraint"""
        self.__weightingScheme = weightingScheme

    @property
    def _experimentalX(self):
        """For internal use only to interface
        ExperimentalConstraint.get_constraints_properties"""
        return self.__experimentalDistances

    @property
    def _experimentalY(self):
        """For internal use only to interface
        ExperimentalConstraint.get_constraints_properties"""
        return self.__experimentalPDF

    @property
    def _modelX(self):
        """For internal use only to interface
        ExperimentalConstraint.get_constraints_properties"""
        return self.__shellCenters

    def listen(self, message, argument=None):
        """
        Listens to any message sent from the Broadcaster.

        :Parameters:
            #. message (object): Any python object to send to constraint's
               listen method.
            #. argument (object): Any type of argument to pass to the listeners.
        """
        if message in ("engine set", "update pdb", "update molecules indexes", "update elements indexes", "update names indexes"):
            if self.engine is not None:
                self.__elementsPairs = sorted(itertools.combinations_with_replacement(self.engine.elements, 2))
            
                # 使用X射线散射因子计算权重
                if self.__weighting == "xray_scattering":
                    self.__weightingScheme = self._calculate_xray_weights()
                    # 为兼容性保留元素权重
                    self._elementsWeight = {}
                    for elem in self.engine.elements:
                        # 使用原子序数作为近似权重（仅用于显示目的）
                        try:
                            self._elementsWeight[elem] = get_element_property.element(elem).atomic_number
                        except:
                            self._elementsWeight[elem] = 1.0
                else:
                    # 原有的权重计算方式
                    self._elementsWeight = get_real_elements_weight(elements=self.engine.elements, 
                                                                   weightsDict=self.__atomsWeight, 
                                                                   weighting=self.__weighting)
                    self.__weightingScheme = get_normalized_weighting(numbers=self.engine.numberOfAtomsPerElement, 
                                                                     weights=self._elementsWeight)
                    for k in self.__weightingScheme:
                        self.__weightingScheme[k] = FLOAT_TYPE(self.__weightingScheme[k])
            else:
                self.__elementsPairs = None
                self._elementsWeight = None
                self.__weightingScheme = None
            
            # 转储到存储库
            self._dump_to_repository({'_PairDistributionConstraint__elementsPairs': self.__elementsPairs,
                                     '_PairDistributionConstraint__weightingScheme': self.__weightingScheme,
                                     '_elementsWeight': self._elementsWeight})
            self.reset_constraint()  # ADDED 2017-JAN-08
        elif message in ("update boundary conditions",):
            self.reset_constraint()

    def set_shape_function_parameters(self, shapeFuncParams):
        """
        Set the shape function. The shape function can be set to None which
        means unsused, or set as a constant shape given by a numpy.ndarray or
        computed from all atoms and updated every 'updateFreq' accepted moves.
        The shape function is subtracted from the total g(r). It must be used
        when non-periodic boundary conditions are used to take into account
        the atomic density drop and to correct for the :math:`\\rho_{0}`
        approximation.

        :Parameters:
            #. shapeFuncParams (None, numpy.ndarray, dict): The shape function
               is  subtracted from the total g(r). It must be used when
               non-periodic boundary conditions are used to take into account
               the atomic density drop and to correct for the :math:`\\rho_{0}`
               approximation. The shape function can be set to None which means
               unsused, or set as a constant shape given by a numpy.ndarray or
               computed from all atoms and updated every 'updateFreq' accepted
               moves. If dict is given the following keywords can be given,
               otherwise default values will be automatically set.\n
               * **rmin (number) default (0.00) :** The minimum distance in
                 :math:`\\AA` considered upon building the histogram prior to
                 computing the shape function. If not defined, rmin will be
                 set automatically to 0.
               * **rmax (None, number) default (None) :** The maximum distance
                 in :math:`\\AA` considered upon building the histogram prior
                 to computing the shape function. If not defnined, rmax will
                 be automatically set to :math:`maximum\ box\ length + 10\\AA`
                 at engine runtime.
               * **dr (number) default (0.5) :** The bin size in :math:`\\AA`
                 considered upon building the histogram prior to computing the
                 shape function. If not defined, it will be automatically set
                 to 0.5.
               * **qmin (number) default (0.001) :** The minimum reciprocal
                 distance q in :math:`\\AA^{-1}` considered to compute the
                 shape function. If not defined, it will be automatically
                 set to 0.001.
               * **qmax (number) default (0.75) :** The maximum reciprocal
                 distance q in :math:`\\AA^{-1}` considered to compute the
                 shape function. If not defined, it will be automatically
                 set to 0.75.
               * **dq (number) default (0.005) :** The reciprocal distance bin
                 size in :math:`\\AA^{-1}` considered to compute the shape
                 function. If not defined, it will be automatically
                 set to 0.005.
               * **updateFreq (integer) default (1000) :** The frequency of
                 recomputing the shape function in number of accpeted moves.
        """
        self._shapeArray = None
        if shapeFuncParams is None:
            self._shapeFuncParams = None
            self._shapeUpdateFreq = 0
        elif isinstance(shapeFuncParams, dict):
            rmin            = FLOAT_TYPE( shapeFuncParams.get('rmin',0.00 ) )
            rmax            =             shapeFuncParams.get('rmax',None )
            dr              = FLOAT_TYPE( shapeFuncParams.get('dr'  ,0.5  ) )
            qmin            = FLOAT_TYPE( shapeFuncParams.get('qmin',0.001) )
            qmax            = FLOAT_TYPE( shapeFuncParams.get('qmax',0.75 ) )
            dq              = FLOAT_TYPE( shapeFuncParams.get('dq'  ,0.005) )
            self._shapeFuncParams = {'rmin':rmin, 'rmax':rmax, 'dr':dr,
                                      'qmin':qmin, 'qmax':qmax, 'dq':dq }
            self._shapeUpdateFreq = INT_TYPE( shapeFuncParams.get('updateFreq',1000) )
        else:
            assert isinstance(shapeFuncParams, (list,tuple,np.ndarray)), LOGGER.error("shapeFuncParams must be None, numpy.ndarray or a dictionary")
            try:
                shapeArray = np.array(shapeFuncParams)
            except:
                raise LOGGER.error("constant shapeFuncParams must be numpy.ndarray castable")
            assert len(shapeFuncParams.shape) == 1, LOGGER.error("numpy.ndarray shapeFuncParams must be of dimension 1")
            assert shapeFuncParams.shape[0] == self.experimentalData.shape[0], LOGGER.error("numpy.ndarray shapeFuncParams must have the same experimental data length")
            for n in shapeFuncParams:
                assert is_number(n), LOGGER.error("numpy.ndarray shapeFuncParams must be numbers")
            self._shapeFuncParams = shapeFuncParams.astype(FLOAT_TYPE)
            self._shapeUpdateFreq = 0
        # dump to repository
        self._dump_to_repository({'_shapeFuncParams': self._shapeFuncParams,
                                  '_shapeUpdateFreq': self._shapeUpdateFreq})

    def set_weighting(self, weighting):
        """
        Set elements weighting. It must be a valid entry of pdbparser atom's
        database.

        :Parameters:
            #. weighting (string): The elements weighting scheme. It must be
               any atomic attribute (atomicNumber, neutronCohb, neutronIncohb,
               neutronCohXs, neutronIncohXs, atomicWeight, covalentRadius)
               defined in pdbparser database. In case of xrays or neutrons
               experimental weights, one can simply set weighting to 'xrays'
               or 'neutrons' and the value will be automatically adjusted to
               respectively 'atomicNumber' and 'neutronCohb'. If attribute
               values are missing in the pdbparser database, atomic weights
               must be given in atomsWeight dictionary argument.
        """
        assert self.engine is None, LOGGER.error("Engine is already set. Reseting weighting is not allowed")
    
        if weighting.lower() in ["xrays","x-rays","xray","x-ray"]:
            LOGGER.fixed("'%s' weighting is set to atomicNumber"%weighting)
            weighting = "atomicNumber"
        elif weighting.lower() in ["neutron","neutrons"]:
            LOGGER.fixed("'%s' weighting is set to neutronCohb"%weighting)
            weighting = "neutronCohb"
        elif weighting.lower() == "xray_scattering":
            LOGGER.fixed("'%s' weighting uses Waasmaier-Kirfel based x-ray scattering factors" % weighting)
            self.__weighting = "xray_scattering"
        else:
            assert is_element_property(weighting), LOGGER.error("weighting is not a valid pdbparser atoms database entry")
            assert weighting != "atomicFormFactor", LOGGER.error("atomicFormFactor weighting is not allowed")
    
        self.__weighting = weighting

    def set_atoms_weight(self, atomsWeight):
        """
        Custom set atoms weight. This is the way to customize setting atoms
        weights different than the given weighting scheme.

        :Parameters:
            #. atomsWeight (None, dict): Atoms weight dictionary where keys are
               atoms element and values are custom weights. If None is given
               or partially given, missing elements weighting will be fully set
               given weighting scheme.
        """
        if atomsWeight is None:
            AW = {}
        else:
            assert isinstance(atomsWeight, dict),LOGGER.error("atomsWeight must be None or a dictionary")
            AW = {}
            for k in atomsWeight:
                assert isinstance(k, basestring),LOGGER.error("atomsWeight keys must be strings")
                try:
                    val = FLOAT_TYPE(atomsWeight[k])
                except:
                    raise LOGGER.error( "atomsWeight values must be numerical")
                AW[k]=val
        # set atomsWeight
        self.__atomsWeight = AW
        # reset weights
        if self.engine is None:
            self.__elementsPairs   = None
            self._elementsWeight   = None
            self.__weightingScheme = None
        else:
            isNormalFrame, isMultiframe, isSubframe = self.engine.get_frame_category(frame=self.engine.usedFrame)
            self.__elementsPairs   = sorted(itertools.combinations_with_replacement(self.engine.elements,2))
            self._elementsWeight   = get_real_elements_weight(elements=self.engine.elements, weightsDict=self.__atomsWeight, weighting=self.__weighting)
            self.__weightingScheme = get_normalized_weighting(numbers=self.engine.numberOfAtomsPerElement, weights=self._elementsWeight)
            for k in self.__weightingScheme:
                self.__weightingScheme[k] = FLOAT_TYPE(self.__weightingScheme[k])
            if isSubframe:
                repo   = self.engine._get_repository()
                assert repo is not None, LOGGER.error("Repository is not defined, not allowed to set atoms weight for a subframe.")
                mframe = self.engine.usedFrame.split(os.sep)[0]
                LOGGER.usage("set_atoms_weight for '%s' subframe. This is going to automatically propagate to all '%s' multiframe subframes."%(self.engine.usedFrame,mframe))
                for subfrm in self.engine.frames[mframe]['frames_name']:
                    frame = os.path.join(mframe,subfrm)
                    if frame != self.engine.usedFrame:
                        elements         = repo.pull(relativePath=os.path.join(frame,'_Engine__elements'))
                        nAtomsPerElement = repo.pull(relativePath=os.path.join(frame,'_Engine__numberOfAtomsPerElement'))
                        elementsWeight   = repo.pull(relativePath=os.path.join(frame,'constraints',self.constraintName,'_elementsWeight'))
                        elementsPairs    = sorted(itertools.combinations_with_replacement(elements,2))
                        elementsWeight   = get_real_elements_weight(elements=elements, weightsDict=self.__atomsWeight, weighting=self.__weighting)
                        weightingScheme  = get_normalized_weighting(numbers=nAtomsPerElement, weights=elementsWeight)
                        for k in self.__weightingScheme:
                            weightingScheme[k] = FLOAT_TYPE(self.__weightingScheme[k])
                    # dump to repository
                    self._dump_to_repository({'_PairDistributionConstraint__elementsPairs'  : elementsPairs,
                                              '_PairDistributionConstraint__weightingScheme': weightingScheme,
                                              '_PairDistributionConstraint__atomsWeight'    : self.__atomsWeight,
                                              '_elementsWeight': elementsWeight},
                                              frame=frame)
            else:
                assert isNormalFrame, LOGGER.error("Not allowed to set_atoms_weight for multiframe")
                # dump to repository
                self._dump_to_repository({'_PairDistributionConstraint__elementsPairs'  : self.__elementsPairs,
                                          '_PairDistributionConstraint__weightingScheme': self.__weightingScheme,
                                          '_PairDistributionConstraint__atomsWeight'    : self.__atomsWeight,
                                          '_elementsWeight': self._elementsWeight})



    def set_window_function(self, windowFunction, frame=None):
        """
        Set convolution window function.

        :Parameters:
             #. windowFunction (None, numpy.ndarray): The window function to
                convolute with the computed pair distribution function of the
                system prior to comparing it with the experimental data.
             #. frame (None, string): Target frame name. If None, engine used
                frame is used. If multiframe is given, all subframes will be
                targeted. If subframe is given, all other multiframe subframes
                will be targeted.
        """
        if windowFunction is not None:
            assert isinstance(windowFunction, np.ndarray), LOGGER.error("windowFunction must be a numpy.ndarray")
            assert windowFunction.dtype.type is FLOAT_TYPE, LOGGER.error("windowFunction type must be %s"%FLOAT_TYPE)
            assert len(windowFunction.shape) == 1, LOGGER.error("windowFunction must be of dimension 1")
            assert len(windowFunction) <= self.experimentalData.shape[0], LOGGER.error("windowFunction length must be smaller than experimental data")
            # normalize window function
            windowFunction /= np.sum(windowFunction)
        # check window size
        # set windowFunction
        self.__windowFunction = windowFunction
        # dump to repository
        usedIncluded, frame, allFrames = get_caller_frames(engine=self.engine,
                                                           frame=frame,
                                                           subframeToAll=True,
                                                           caller="%s.%s"%(self.__class__.__name__,inspect.stack()[0][3]) )
        if usedIncluded:
            self.__windowFunction = windowFunction
        for frm in allFrames:
            self._dump_to_repository({'_PairDistributionConstraint__windowFunction': self.__windowFunction}, frame=frm)

    def set_experimental_data(self, experimentalData):
        """
        Set constraint's experimental data.

        :Parameters:
            #. experimentalData (numpy.ndarray, string): The experimental
               data as numpy.ndarray or string path to load data using
               numpy.loadtxt function.
        """
        # get experimental data
        super(PairDistributionConstraint, self).set_experimental_data(experimentalData=experimentalData)
        self.__bin = FLOAT_TYPE(self.experimentalData[1,0] - self.experimentalData[0,0])
        # dump to repository
        self._dump_to_repository({'_PairDistributionConstraint__bin': self.__bin})
        # set limits
        self.set_limits(self.limits)


    def set_limits(self, limits):
        """
        Set the histogram computation limits.

        :Parameters:
            #. limits (None, tuple, list): Distance limits to bound
               experimental data and compute histograms.
               If None is given, the limits will be automatically
               set the the min and max distance of the experimental data.
               Otherwise, a tuple of exactly two items where the first is the
               minimum distance or None and the second is the maximum distance
               or None.
        """
        self._ExperimentalConstraint__set_limits(limits)
        # set minimumDistance, maximumDistance
        self.__minimumDistance = FLOAT_TYPE(self.experimentalData[self.limitsIndexStart,0] - self.__bin/2. )
        self.__maximumDistance = FLOAT_TYPE(self.experimentalData[self.limitsIndexEnd ,0]  + self.__bin/2. )
        self.__shellCenters    = np.array([self.experimentalData[idx,0] for idx in range(self.limitsIndexStart,self.limitsIndexEnd +1)],dtype=FLOAT_TYPE)
        # set histogram edges
        edges = [self.experimentalData[idx,0] - self.__bin/2. for idx in range(self.limitsIndexStart,self.limitsIndexEnd +1)]
        edges.append( self.experimentalData[self.limitsIndexEnd ,0] + self.__bin/2. )
        self.__edges = np.array(edges, dtype=FLOAT_TYPE)
        # set histogram size
        self.__histogramSize = INT_TYPE( len(self.__edges)-1 )
        # set shell centers and volumes
        self.__shellVolumes = FLOAT_TYPE(4.0/3.)*PI*((self.__edges[1:])**3 - self.__edges[0:-1]**3)
        # set experimental distances and pdf
        self.__experimentalDistances = self.experimentalData[self.limitsIndexStart:self.limitsIndexEnd +1,0]
        self.__experimentalPDF       = self.experimentalData[self.limitsIndexStart:self.limitsIndexEnd +1,1]
        # dump to repository
        self._dump_to_repository({'_PairDistributionConstraint__minimumDistance'      : self.__minimumDistance,
                                  '_PairDistributionConstraint__maximumDistance'      : self.__maximumDistance,
                                  '_PairDistributionConstraint__shellCenters'         : self.__shellCenters,
                                  '_PairDistributionConstraint__edges'                : self.__edges,
                                  '_PairDistributionConstraint__histogramSize'        : self.__histogramSize,
                                  '_PairDistributionConstraint__shellVolumes'         : self.__shellVolumes,
                                  '_PairDistributionConstraint__experimentalDistances': self.__experimentalDistances,
                                  '_PairDistributionConstraint__experimentalPDF'      : self.__experimentalPDF})
        # set used dataWeights
        self._set_used_data_weights(limitsIndexStart=self.limitsIndexStart, limitsIndexEnd=self.limitsIndexEnd )
        # reset constraint
        self.reset_constraint()

    def check_experimental_data(self, experimentalData):
        """
        Check whether experimental data is correct.

        :Parameters:
            #. experimentalData (object): Experimental data to check.

        :Returns:
            #. result (boolean): Whether it is correct or not.
            #. message (str): Checking message that explains whats's
               wrong with the given data.
        """
        if not isinstance(experimentalData, np.ndarray):
            return False, "experimentalData must be a numpy.ndarray"
        if experimentalData.dtype.type is not FLOAT_TYPE:
            return False, "experimentalData type must be %s"%FLOAT_TYPE
        if len(experimentalData.shape) !=2:
            return False, "experimentalData must be of dimension 2"
        if experimentalData.shape[1] !=2:
            return False, "experimentalData must have only 2 columns"
        # check distances order
        inOrder = (np.array(sorted(experimentalData[:,0]), dtype=FLOAT_TYPE)-experimentalData[:,0])<=PRECISION
        if not np.all(inOrder):
            return False, "experimentalData distances are not sorted in order"
        if experimentalData[0][0]<0:
            return False, "experimentalData distances min value is found negative"
        bin  = experimentalData[1,0] -experimentalData[0,0]
        bins = experimentalData[1:,0]-experimentalData[0:-1,0]
        for b in bins:
            if np.abs(b-bin)>PRECISION:
                return False, "experimentalData distances bins are found not coherent"
        # data format is correct
        return True, ""

    def compute_standard_error(self, modelData):
        """
        Compute the standard error (StdErr) as the squared deviations
        between model computed data and the experimental ones.

        .. math::
            StdErr = \\sum \\limits_{i}^{N} W_{i}(Y(X_{i})-F(X_{i}))^{2}

        Where:\n
        :math:`N` is the total number of experimental data points. \n
        :math:`W_{i}` is the data point weight. It becomes equivalent to 1 when dataWeights is set to None. \n
        :math:`Y(X_{i})` is the experimental data point :math:`X_{i}`. \n
        :math:`F(X_{i})` is the computed from the model data  :math:`X_{i}`. \n

        :Parameters:
            #. modelData (numpy.ndarray): The data to compare with the
               experimental one and compute the squared deviation.

        :Returns:
            #. standardError (number): The calculated constraint's
               standardError.
        """
        diff = self.__experimentalPDF-modelData
        # return squared deviation
        if self._usedDataWeights is None:
            return np.add.reduce((diff)**2)
        else:
            return np.add.reduce(self._usedDataWeights*((diff)**2))

    def update_standard_error(self):
        """ Compute and set constraint's standardError."""
        # set standardError
        totalPDF = self.get_constraint_value()["total"]
        self.set_standard_error(self.compute_standard_error(modelData = totalPDF))


    def __get_total_Gr(self, data, rho0=None):
        """
        计算g(r)并直接用于约束。保持方法名不变以兼容现有代码。
        
        :Parameters:
            #. data (dict): 含有intra和inter直方图数据的字典
            #. rho0 (float): 现在这个参数不再使用，但为了与调用代码兼容保留它
            
        :Returns:
            #. gr (numpy.ndarray): 计算的g(r)数组
        """
        # rho0参数在此版本中不再使用，但为了保持兼容性仍然接受它
        
        # 初始化 gr 数组
        gr = np.zeros(self.__histogramSize, dtype=FLOAT_TYPE)
        for pair in self.__elementsPairs:
            # 获取权重方案
            wij = self.__weightingScheme.get(pair[0]+"-"+pair[1], None)
            if wij is None:
                wij = self.__weightingScheme[pair[1]+"-"+pair[0]]
            # 获取每种元素的原子数
            ni = self.engine.numberOfAtomsPerElement[pair[0]]
            nj = self.engine.numberOfAtomsPerElement[pair[1]]
            # 获取元素索引
            idi = self.engine.elements.index(pair[0])
            idj = self.engine.elements.index(pair[1])
            # 获取 Nij
            if idi == idj:
                Nij = ni*(ni-1)/2.0
                Dij = FLOAT_TYPE( Nij/self.engine.volume )
                nij = data["intra"][idi,idj,:]+data["inter"][idi,idj,:]
                gr += wij*nij/Dij
            else:
                Nij = ni*nj
                Dij = FLOAT_TYPE( Nij/self.engine.volume )
                nij = data["intra"][idi,idj,:]+data["intra"][idj,idi,:] + data["inter"][idi,idj,:]+data["inter"][idj,idi,:]
                gr += wij*nij/Dij
        # 除以壳层体积
        gr /= self.shellVolumes
        
        # 注意：不再将g(r)转换为G(r)
        # 之前的代码: Gr = (4.*PI*self.__shellCenters*rho0)*(Gr-1)
        
        # 移除形状函数
        if self._shapeArray is not None:
            gr -= self._shapeArray
        
        # 乘以缩放因子
        self._fittedScaleFactor = self.get_adjusted_scale_factor(self.experimentalPDF, gr, self._usedDataWeights)
        if self._fittedScaleFactor != 1:
            gr *= FLOAT_TYPE(self._fittedScaleFactor)
        
        # 应用多框架先验和权重
        gr = self._apply_multiframe_prior(gr)
        
        # 与窗口函数卷积
        if self.__windowFunction is not None:
            gr = np.convolve(gr, self.__windowFunction, 'same')
        
        # 返回数组
        return gr


    def _get_constraint_value(self, data, applyMultiframePrior=True):
        """
        计算约束数据（intra, inter, total, total_no_window）
        缩放因子会被应用，但多框架先验不会

        :Parameters:
            #. data (dict): 计算约束值所需的数据字典
            #. applyMultiframePrior (boolean): 是否应用子框架的权重和先验
            
        :Returns:
            #. output (dict): 包含计算结果的字典
        """
        output = {}
        for pair in self.__elementsPairs:
            output["rdf_intra_%s-%s" % pair] = np.zeros(self.__histogramSize, dtype=FLOAT_TYPE)
            output["rdf_inter_%s-%s" % pair] = np.zeros(self.__histogramSize, dtype=FLOAT_TYPE)
            output["rdf_total_%s-%s" % pair] = np.zeros(self.__histogramSize, dtype=FLOAT_TYPE)
        gr = np.zeros(self.__histogramSize, dtype=FLOAT_TYPE)
        for pair in self.__elementsPairs:
            # get weighting scheme
            wij = self.__weightingScheme.get(pair[0]+"-"+pair[1], None)
            if wij is None:
                wij = self.__weightingScheme[pair[1]+"-"+pair[0]]
            # get number of atoms per element
            ni = self.engine.numberOfAtomsPerElement[pair[0]]
            nj = self.engine.numberOfAtomsPerElement[pair[1]]
            # get index of element
            idi = self.engine.elements.index(pair[0])
            idj = self.engine.elements.index(pair[1])
            # get Nij
            if idi == idj:
                Nij = FLOAT_TYPE( ni*(ni-1)/2.0 )
                output["rdf_intra_%s-%s" % pair] += data["intra"][idi,idj,:]
                output["rdf_inter_%s-%s" % pair] += data["inter"][idi,idj,:]
            else:
                Nij = FLOAT_TYPE( ni*nj )
                output["rdf_intra_%s-%s" % pair] += data["intra"][idi,idj,:] + data["intra"][idj,idi,:]
                output["rdf_inter_%s-%s" % pair] += data["inter"][idi,idj,:] + data["inter"][idj,idi,:]
            # compute g(r)
            nij = output["rdf_intra_%s-%s" % pair] + output["rdf_inter_%s-%s" % pair]
            dij = nij/self.__shellVolumes
            Dij = Nij/self.engine.volume
            gr += wij*dij/Dij
            # calculate intensityFactor
            intensityFactor = (self.engine.volume*wij)/(Nij*self.__shellVolumes)
            # divide by factor
            output["rdf_intra_%s-%s" % pair] *= intensityFactor
            output["rdf_inter_%s-%s" % pair] *= intensityFactor
            output["rdf_total_%s-%s" % pair]  = output["rdf_intra_%s-%s" % pair] + output["rdf_inter_%s-%s" % pair]

        # 存储g(r)而非G(r)
        output["total_no_window"] = gr

        # 移除形状函数
        if self._shapeArray is not None:
            output["total_no_window"] -= self._shapeArray
            
        # multiply by scale factor
        if self.scaleFactor != 1:
            output["total_no_window"] *= self.scaleFactor
            
        # apply multiframe prior and weight
        if applyMultiframePrior:
            output["total_no_window"] = self._apply_multiframe_prior(output["total_no_window"])
            
        # convolve total with window function
        if self.__windowFunction is not None:
            output["total"] = np.convolve(output["total_no_window"], self.__windowFunction, 'same')
        else:
            output["total"] = output["total_no_window"]
            
        return output

    def get_constraint_value(self, applyMultiframePrior=True):
        """
        计算所有偏对分布函数(PDFs)。

        :Parameters:
            #. applyMultiframePrior (boolean): 是否将子框架权重和先验应用于总量。
               只有当使用的框架是子框架且定义了子框架权重和先验时才有效。

        :Returns:
            #. PDFs (dictionary): PDFs字典，键是元素的分子内和分子间PDFs，
               值是计算的PDFs。
        """
        if self.data is None:
            LOGGER.warn("data must be computed first using 'compute_data' method.")
            return {}
        return self._get_constraint_value(self.data, applyMultiframePrior=applyMultiframePrior)

    def get_constraint_original_value(self):
        """
        计算所有偏对分布函数(PDFs)。

        :Returns:
            #. PDFs (dictionary): PDFs字典，键是元素的分子内和分子间PDFs，
               值是计算的PDFs。
        """
        if self.originalData is None:
            LOGGER.warn("originalData must be computed first using 'compute_data' method.")
            return {}
        return self._get_constraint_value(self.originalData)

    @reset_if_collected_out_of_date
    def compute_data(self, update=True):
        """ 
        计算约束数据。

        :Parameters:
            #. update (boolean): 是否用新计算更新约束数据和标准误差。
               如果在随机引擎运行时数据由其他线程或进程计算和更新，可能会
               导致约束状态改变，从而导致运行中不再接受额外的移动。

        :Returns:
            #. data (dict): 约束数据字典
            #. standardError (float): 约束标准误差
        """
        intra,inter = full_pairs_histograms_coords( boxCoords        = self.engine.boxCoordinates,
                                                    basis            = self.engine.basisVectors,
                                                    isPBC            = self.engine.isPBC,
                                                    moleculeIndex    = self.engine.moleculesIndex,
                                                    elementIndex     = self.engine.elementsIndex,
                                                    numberOfElements = self.engine.numberOfElements,
                                                    minDistance      = self.__minimumDistance,
                                                    maxDistance      = self.__maximumDistance,
                                                    histSize         = self.__histogramSize,
                                                    bin              = self.__bin,
                                                    ncores           = self.engine._runtime_ncores )
        # create data and compute standard error
        data     = {"intra":intra, "inter":inter}
        totalPDF = self.__get_total_Gr(data, rho0=self.engine.numberDensity)
        stdError = self.compute_standard_error(modelData = totalPDF)
        # update
        if update:
            self.set_data(data)
            self.set_active_atoms_data_before_move(None)
            self.set_active_atoms_data_after_move(None)
            # set standardError
            self.set_standard_error(stdError)
            # set original data
            if self.originalData is None:
                self._set_original_data(self.data)
        # return
        return data, stdError

    def compute_before_move(self, realIndexes, relativeIndexes):
        """
        在执行移动前计算约束。

        :Parameters:
            #. realIndexes (numpy.ndarray): 此处不使用。
            #. relativeIndexes (numpy.ndarray): 将应用移动的组原子相对索引。
        """
        intraM,interM = multiple_pairs_histograms_coords( indexes          = relativeIndexes,
                                                          boxCoords        = self.engine.boxCoordinates,
                                                          basis            = self.engine.basisVectors,
                                                          isPBC            = self.engine.isPBC,
                                                          moleculeIndex    = self.engine.moleculesIndex,
                                                          elementIndex     = self.engine.elementsIndex,
                                                          numberOfElements = self.engine.numberOfElements,
                                                          minDistance      = self.__minimumDistance,
                                                          maxDistance      = self.__maximumDistance,
                                                          histSize         = self.__histogramSize,
                                                          bin              = self.__bin,
                                                          allAtoms         = True,
                                                          ncores           = self.engine._runtime_ncores)
        intraF,interF = full_pairs_histograms_coords( boxCoords        = self.engine.boxCoordinates[relativeIndexes],
                                                      basis            = self.engine.basisVectors,
                                                      isPBC            = self.engine.isPBC,
                                                      moleculeIndex    = self.engine.moleculesIndex[relativeIndexes],
                                                      elementIndex     = self.engine.elementsIndex[relativeIndexes],
                                                      numberOfElements = self.engine.numberOfElements,
                                                      minDistance      = self.__minimumDistance,
                                                      maxDistance      = self.__maximumDistance,
                                                      histSize         = self.__histogramSize,
                                                      bin              = self.__bin,
                                                      ncores           = self.engine._runtime_ncores )
        # set active atoms data
        self.set_active_atoms_data_before_move( {"intra":intraM-intraF, "inter":interM-interF} )
        self.set_active_atoms_data_after_move(None)

    def compute_after_move(self, realIndexes, relativeIndexes, movedBoxCoordinates):
        """
        在执行移动后计算约束

        :Parameters:
            #. realIndexes (numpy.ndarray): 此处不使用。
            #. relativeIndexes (numpy.ndarray): 将应用移动的组原子相对索引。
            #. movedBoxCoordinates (numpy.ndarray): 移动后原子的新坐标。
        """
        # change coordinates temporarily
        boxData = np.array(self.engine.boxCoordinates[relativeIndexes], dtype=FLOAT_TYPE)
        self.engine.boxCoordinates[relativeIndexes] = movedBoxCoordinates
        # calculate pair distribution function
        intraM,interM = multiple_pairs_histograms_coords( indexes          = relativeIndexes,
                                                          boxCoords        = self.engine.boxCoordinates,
                                                          basis            = self.engine.basisVectors,
                                                          isPBC            = self.engine.isPBC,
                                                          moleculeIndex    = self.engine.moleculesIndex,
                                                          elementIndex     = self.engine.elementsIndex,
                                                          numberOfElements = self.engine.numberOfElements,
                                                          minDistance      = self.__minimumDistance,
                                                          maxDistance      = self.__maximumDistance,
                                                          histSize         = self.__histogramSize,
                                                          bin              = self.__bin,
                                                          allAtoms         = True,
                                                          ncores           = self.engine._runtime_ncores )
        intraF,interF = full_pairs_histograms_coords( boxCoords        = self.engine.boxCoordinates[relativeIndexes],
                                                      basis            = self.engine.basisVectors,
                                                      isPBC            = self.engine.isPBC,
                                                      moleculeIndex    = self.engine.moleculesIndex[relativeIndexes],
                                                      elementIndex     = self.engine.elementsIndex[relativeIndexes],
                                                      numberOfElements = self.engine.numberOfElements,
                                                      minDistance      = self.__minimumDistance,
                                                      maxDistance      = self.__maximumDistance,
                                                      histSize         = self.__histogramSize,
                                                      bin              = self.__bin,
                                                      ncores           = self.engine._runtime_ncores )
        # set active atoms data
        self.set_active_atoms_data_after_move( {"intra":intraM-intraF, "inter":interM-interF} )
        # reset coordinates
        self.engine.boxCoordinates[relativeIndexes] = boxData
        # compute and set standardError after move
        dataIntra = self.data["intra"]-self.activeAtomsDataBeforeMove["intra"]+self.activeAtomsDataAfterMove["intra"]
        dataInter = self.data["inter"]-self.activeAtomsDataBeforeMove["inter"]+self.activeAtomsDataAfterMove["inter"]
        totalPDF  = self.__get_total_Gr({"intra":dataIntra, "inter":dataInter}, rho0=self.engine.numberDensity)
        self.set_after_move_standard_error( self.compute_standard_error(modelData = totalPDF) )
        # increment tried
        self.increment_tried()

    def accept_move(self, realIndexes, relativeIndexes):
        """
        接受移动

        :Parameters:
            #. realIndexes (numpy.ndarray): 此处不使用。
            #. relativeIndexes (numpy.ndarray): 此处不使用。
        """
        dataIntra = self.data["intra"]-self.activeAtomsDataBeforeMove["intra"]+self.activeAtomsDataAfterMove["intra"]
        dataInter = self.data["inter"]-self.activeAtomsDataBeforeMove["inter"]+self.activeAtomsDataAfterMove["inter"]
        # change permanently _data
        self.set_data( {"intra":dataIntra, "inter":dataInter} )
        # reset activeAtoms data
        self.set_active_atoms_data_before_move(None)
        self.set_active_atoms_data_after_move(None)
        # update standardError
        self.set_standard_error( self.afterMoveStandardError )
        self.set_after_move_standard_error( None )
        # set new scale factor
        self._set_fitted_scale_factor_value(self._fittedScaleFactor)
        # increment accepted
        self.increment_accepted()

    def reject_move(self, realIndexes, relativeIndexes):
        """
        拒绝移动

        :Parameters:
            #. realIndexes (numpy.ndarray): 此处不使用。
            #. relativeIndexes (numpy.ndarray): 此处不使用。
        """
        # reset activeAtoms data
        self.set_active_atoms_data_before_move(None)
        self.set_active_atoms_data_after_move(None)
        # update standardError
        self.set_after_move_standard_error( None )

    def compute_as_if_amputated(self, realIndex, relativeIndex):
        """
        计算并返回约束的数据和标准误差，如同给定原子被截除一样。

        :Parameters:
            #. realIndex (numpy.ndarray): 作为单一元素的numpy数组的原子索引。
            #. relativeIndex (numpy.ndarray): 作为单一元素的numpy数组的原子相对索引。
        """
        # compute data
        self.compute_before_move(realIndexes=realIndex, relativeIndexes=relativeIndex)
        dataIntra = self.data["intra"]-self.activeAtomsDataBeforeMove["intra"]
        dataInter = self.data["inter"]-self.activeAtomsDataBeforeMove["inter"]
        data      = {"intra":dataIntra, "inter":dataInter}
        # temporarily adjust self.__weightingScheme
        weightingScheme = self.__weightingScheme
        relativeIndex = relativeIndex[0]
        selectedElement = self.engine.allElements[relativeIndex]
        self.engine.numberOfAtomsPerElement[selectedElement] -= 1
        self.__weightingScheme = get_normalized_weighting(numbers=self.engine.numberOfAtomsPerElement, weights=self._elementsWeight )
        for k in self.__weightingScheme:
            self.__weightingScheme[k] = FLOAT_TYPE(self.__weightingScheme[k])
        # compute standard error
        if not self.engine._RT_moveGenerator.allowFittingScaleFactor:
            SF = self.adjustScaleFactorFrequency
            self._set_adjust_scale_factor_frequency(0)
        
        # 计算原子移除后的密度（仅用于传参，现在不再使用）
        rho0 = ((self.engine.numberOfAtoms-1)/self.engine.volume).astype(FLOAT_TYPE)
        totalPDF = self.__get_total_Gr(data, rho0=rho0)
        
        standardError = self.compute_standard_error(modelData = totalPDF)
        if not self.engine._RT_moveGenerator.allowFittingScaleFactor:
            self._set_adjust_scale_factor_frequency(SF)
        # reset activeAtoms data
        self.set_active_atoms_data_before_move(None)
        # set data
        self.set_amputation_data( {'data':data, 'weightingScheme':self.__weightingScheme} )
        # compute standard error
        self.set_amputation_standard_error( standardError )
        # reset weightingScheme and number of atoms per element
        self.__weightingScheme = weightingScheme
        self.engine.numberOfAtomsPerElement[selectedElement] += 1
        print(self.engine.numberOfAtoms, rho0, self.engine.numberOfAtomsPerElement, self.engine.numberDensity)

    def accept_amputation(self, realIndex, relativeIndex):
        """
        接受截除原子，并相应地设置约束数据和标准误差。

        :Parameters:
            #. realIndex (numpy.ndarray): 此处不使用。
            #. relativeIndex (numpy.ndarray): 此处不使用。
        """
        self.set_data( self.amputationData['data'] )
        self.__weightingScheme = self.amputationData['weightingScheme']
        self.set_standard_error( self.amputationStandardError )
        self.set_amputation_data( None )
        self.set_amputation_standard_error( None )
        # set new scale factor
        self._set_fitted_scale_factor_value(self._fittedScaleFactor)

    def reject_amputation(self, realIndex, relativeIndex):
        """
        拒绝截除原子，并相应地设置约束数据和标准误差。

        :Parameters:
            #. realIndex (numpy.ndarray): 此处不使用。
            #. relativeIndex (numpy.ndarray): 此处不使用。
        """
        self.set_amputation_data( None )
        self.set_amputation_standard_error( None )

    def _on_collector_collect_atom(self, realIndex):
        pass

    def _on_collector_release_atom(self, realIndex):
        pass

    def get_multiframe_weights(self, frame):
        """
        获取多框架权重。

        :Parameters:
            #. frame (str): 框架名称。

        :Returns:
            #. weights (OrderedDict): 包含子框架权重的有序字典。
        """
        from collections import OrderedDict
        isNormalFrame, isMultiframe, isSubframe = self.engine.get_frame_category(frame=frame)
        assert isMultiframe, LOGGER.error("Given frame '%s' is not a multiframe"%frame)
        repo    = self.engine._get_repository()
        weights = OrderedDict()
        for frm in self.engine.frames[frame]['frames_name']:
            value = repo.pull(relativePath=os.path.join(frame,frm,'constraints',self.constraintName,'_ExperimentalConstraint__multiframeWeight'))
            weights[frm] = value
        return weights

    def _constraint_copy_needs_lut(self):
        return {'_PairDistributionConstraint__elementsPairs'        :'_PairDistributionConstraint__elementsPairs',
                '_PairDistributionConstraint__histogramSize'        :'_PairDistributionConstraint__histogramSize',
                '_PairDistributionConstraint__weightingScheme'      :'_PairDistributionConstraint__weightingScheme',
                '_PairDistributionConstraint__shellVolumes'         :'_PairDistributionConstraint__shellVolumes',
                '_PairDistributionConstraint__shellCenters'         :'_PairDistributionConstraint__shellCenters',
                '_PairDistributionConstraint__windowFunction'       :'_PairDistributionConstraint__windowFunction',
                '_PairDistributionConstraint__experimentalDistances':'_PairDistributionConstraint__experimentalDistances',
                '_PairDistributionConstraint__experimentalPDF'      :'_PairDistributionConstraint__experimentalPDF',
                '_PairDistributionConstraint__minimumDistance'      :'_PairDistributionConstraint__minimumDistance',
                '_PairDistributionConstraint__maximumDistance'      :'_PairDistributionConstraint__maximumDistance',
                '_PairDistributionConstraint__bin'                  :'_PairDistributionConstraint__bin',
                '_q_range'                                          :'_q_range',
                '_shapeArray'                                       :'_shapeArray',
                '_ExperimentalConstraint__scaleFactor'              :'_ExperimentalConstraint__scaleFactor',
                '_ExperimentalConstraint__dataWeights'              :'_ExperimentalConstraint__dataWeights',
                '_ExperimentalConstraint__multiframePrior'          :'_ExperimentalConstraint__multiframePrior',
                '_ExperimentalConstraint__multiframeWeight'         :'_ExperimentalConstraint__multiframeWeight',
                '_ExperimentalConstraint__limits'                   :'_ExperimentalConstraint__limits',
                '_ExperimentalConstraint__limitsIndexStart'         :'_ExperimentalConstraint__limitsIndexStart',
                '_ExperimentalConstraint__limitsIndexEnd'           :'_ExperimentalConstraint__limitsIndexEnd',
                '_usedDataWeights'                                  :'_usedDataWeights',
                '_Constraint__used'                                 :'_Constraint__used',
                '_Constraint__data'                                 :'_Constraint__data',
                '_Constraint__state'                                :'_Constraint__state',
                '_Engine__state'                                    :'_Engine__state',
                '_Engine__boxCoordinates'                           :'_Engine__boxCoordinates',
                '_Engine__basisVectors'                             :'_Engine__basisVectors',
                '_Engine__isPBC'                                    :'_Engine__isPBC',
                '_Engine__moleculesIndex'                           :'_Engine__moleculesIndex',
                '_Engine__elementsIndex'                            :'_Engine__elementsIndex',
                '_Engine__numberOfAtomsPerElement'                  :'_Engine__numberOfAtomsPerElement',
                '_Engine__elements'                                 :'_Engine__elements',
                '_Engine__numberDensity'                            :'_Engine__numberDensity',
                '_Engine__volume'                                   :'_Engine__volume',
                '_atomsCollector'                                   :'_atomsCollector',
                ('engine','_atomsCollector')                        :'_atomsCollector',
               }

    def plot(self, xlabelParams={'xlabel':'$r(\\AA)$', 'size':10},
                   ylabelParams={'ylabel':'$g(r)$', 'size':10},
                   **kwargs):
        """
        ExperimentalConstraint.plot的别名，带有额外参数

        :Additional/Adjusted Parameters:
            #. xlabelParams (None, dict): 修改的matplotlib.axes.Axes.set_xlabel参数
            #. ylabelParams (None, dict): 修改的matplotlib.axes.Axes.set_ylabel参数
        """
        return super(PairDistributionConstraint, self).plot(xlabelParams=xlabelParams,
                                                           ylabelParams=ylabelParams,
                                                           **kwargs)
