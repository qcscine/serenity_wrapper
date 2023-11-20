__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import scine_utilities as utils

manager = utils.core.ModuleManager.get_instance()
current_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.dirname(os.path.dirname(os.path.dirname(current_path)))
os.environ['SERENITY_RESOURCES'] = os.path.join(current_path, 'data/')
if not manager.module_loaded('Serenity'):
    if os.path.exists(os.path.join(current_path, 'serenity.module.so')):
        manager.load(os.path.join(current_path, 'serenity.module.so'))
        if not manager.module_loaded('Serenity'):
            raise ImportError('The serenity.module.so was found but could not be loaded.')
    elif os.path.exists(os.path.join(lib_path, 'serenity.module.so')):
        manager.load(os.path.join(lib_path, 'serenity.module.so'))
        if not manager.module_loaded('Serenity'):
            raise ImportError('The serenity.module.so was found but could not be loaded.')
    else:
        raise ImportError('The serenity.module.so could not be located.')
