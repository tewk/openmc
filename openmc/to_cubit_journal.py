from openmc.region import Region, Complement, Intersection, Union
from openmc.surface import Halfspace
from openmc.lattice import Lattice

def to_cubit_journal(geom, seen=set(), world=[60,60,60], filename=None, to_cubit=False ):
    w = world
    cid = 1
    def lastid():
        nonlocal cid 
        id = cid
        cid = cid + 1
        return id
    cmds = []

    def make_world_brick():
        pass
        
    def surface_to_cubit_journal(node, w, indent = 0):
        def ind():
            return ' ' * (2*indent)
        if isinstance(node, Halfspace):
            #if not ( node.surface in seen ):
                seen.add( node.surface )
                surface = node.surface
               # print( ind(), "Surface:", surface )

                def reverse():
                    if surface.name.count( "maximum" ) > 0 : #or node.side != '-'
                        return "reverse"
                    else:
                        return ""

                if surface._type == "plane":
                    raise "plane not implemented"
                    pass
                elif surface._type == "x-plane":
                    cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                    id = lastid()
                    cmds.append( f"section volume {id} with xplane offset {surface.coefficients['x0']} {reverse()}")
                    return id
                elif surface._type == "y-plane":
                    cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                    id = lastid()
                    cmds.append( f"section volume {id} with yplane offset {surface.coefficients['y0']} {reverse()}")
                    return id
                elif surface._type == "z-plane":
                    cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                    id = lastid()
                    cmds.append( f"section volume {id} with zplane offset {surface.coefficients['z0']} {reverse()}")
                    return id
                elif surface._type == "cylinder":
                    raise "cylinder not implemented"
                elif surface._type == "x-cylinder":
                    cmds.append( f"cylinder height {w[0]} radius {surface.coefficients['r']}")
                    id = lastid()
                    cmds.append( f"rotate volume {id} about y angle 90")
                    if node.side != '-':
                        cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                        wid = lastid()
                        cmds.append( f"subtract vol {id} from vol {wid}" )
                        return wid
                    return id
                elif surface._type == "y-cylinder":
                    cmds.append( f"cylinder height {w[1]} radius {surface.coefficients['r']}")
                    id = lastid()
                    cmds.append( f"rotate volume {id} about x angle 90")
                    if node.side != '-':
                        cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                        wid = lastid()
                        cmds.append( f"subtract vol {id} from vol {wid}" )
                        return wid
                    return id
                elif surface._type == "z-cylinder":
                    cmds.append( f"cylinder height {w[2]} radius {surface.coefficients['r']}")
                    id = lastid()
                    if node.side != '-':
                        cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                        wid = lastid()
                        cmds.append( f"subtract vol {id} from vol {wid}" )
                        return wid
                    return id
                elif surface._type == "sphere":
                    print( ind(), f"surface radius {surface._radius}" )
                    pass
                elif surface._type == "cone":
                    raise "cone not implemented"
                    pass
                elif surface._type == "x-cone":
                    cmds.append( f"create frustum height {w[2]} radius {surface.coefficients['r']} top {surface.coefficients['r']}")
                    id = lastid()
                    cmds.append( f"rotate volume {id} about y angle 90")
                    if node.side != '-':
                        cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                        wid = lastid()
                        cmds.append( f"subtract vol {id} from vol {wid}" )
                        return wid
                    return id
                elif surface._type == "y-cone":
                    cmds.append( f"create frustum height {w[2]} radius {surface.coefficients['r']} top {surface.coefficients['r']}")
                    id = lastid()
                    cmds.append( f"rotate volume {id} about x angle 90")
                    if node.side != '-':
                        cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                        wid = lastid()
                        cmds.append( f"subtract vol {id} from vol {wid}" )
                        return wid
                    return id
                elif surface._type == "z-cone":
                    cmds.append( f"create frustum height {w[2]} radius {surface.coefficients['r']} top {surface.coefficients['r']}")
                    id = lastid()
                    if node.side != '-':
                        cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                        wid = lastid()
                        cmds.append( f"subtract vol {id} from vol {wid}" )
                        return wid
                    return id
                elif surface._type == "x-torus":
                    cmds.append( f"torus major radius {surface.coefficients['r']} minor radius {surface.coefficients['r']}")
                    id = lastid()
                    cmds.append( f"rotate volume {id} about y angle 90")
                    if node.side != '-':
                        cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                        wid = lastid()
                        cmds.append( f"subtract vol {id} from vol {wid}" )
                        return wid
                    return id
                elif surface._type == "y-torus":
                    cmds.append( f"torus major radius {surface.coefficients['r']} minor radius {surface.coefficients['r']}")
                    id = lastid()
                    cmds.append( f"rotate volume {id} about x angle 90")
                    if node.side != '-':
                        cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                        wid = lastid()
                        cmds.append( f"subtract vol {id} from vol {wid}" )
                        return wid
                    return id
                elif surface._type == "z-torus":
                    cmds.append( f"torus major radius {surface.coefficients['r']} minor radius {surface.coefficients['r']}")
                    id = lastid()
                    if node.side != '-':
                        cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                        wid = lastid()
                        cmds.append( f"subtract vol {id} from vol {wid}" )
                        return wid
                    return id
                else:
                    raise f"{surface.type} not implemented"
            #else:
            #    print( ind(), node )

        elif isinstance(node, Complement):
            #print( ind(), "Complement:" )
            return surface_to_cubit_journal(node.node, w, indent + 1)
        elif isinstance(node, Intersection):
            #print( ind(), "Intersection:" )
            surfaces = []
            for subnode in node:
                #print( ind(), "Subnode:", subnode )
                id = surface_to_cubit_journal( subnode, w, indent + 1)
                surfaces.append( id )
            if len( surfaces) > 1:
                last = surfaces[0]
                for s in surfaces[1:]:
                    cmds.append( f"intersect {last} {s}" )
                    last = s
            return surfaces[-1]
        elif isinstance(node, Union):
            #print( ind(), "Union:" )
            surfaces = []
            for subnode in node:
                ##print( ind(), "Subnode:", subnode )
                id = surface_to_cubit_journal( subnode, w, indent + 1)
                surfaces.append( id )
            if len( surfaces) > 1:
                last = surfaces[0]
                for s in surfaces[1:]:
                    cmds.append( f"unite {last} {s}" )
                    last = s
            return surfaces[-1]
        else:
            ##print( ind(), "Composite:", node )
            for subnode in node:
                ##print( ind(), "Subnode:", subnode )
                surface_to_cubit_journal( subnode, w, indent + 1)

    def process_node_or_fill( node, w, indent = 0 ):
        def ind():
            return ' ' * (2*indent)
        #if not ( node in seen ):
        #print(  ind(), "Node:", node )
        seen.add( node )
        results = []
        if hasattr( node, "region" ):
            #print(  ind(), "Region:", node.region )
            id = surface_to_cubit_journal( node.region, w, indent )
            results.append( id )
        if hasattr( node, "fill" ) and isinstance(node.fill, Lattice):
            #print(  ind(), "Fill:", node.fill )
            id = process_node_or_fill( node.fill, w, indent + 1 )
            results.append( id )
        if hasattr( node, "universes" ):
            #print(  ind(), "Universes:", node )
            pitch = node._pitch
            ll = [ node.lower_left[0], node.lower_left[1] ]
            ll[0] = ll[0] + pitch[0] / 2
            ll[1] = ll[1] + pitch[1] / 2
            i = 0
            for us in node.universes:
                j = 0
                for u in us:
                    for n, cell in u._cells.items():
                        #print(  ind(), "UCell:", n, cell )
                        id = process_node_or_fill( cell, [pitch[0], pitch[1], w[2] ], indent + 1 )
                        ids = str( id )
                        if isinstance( id, list ):
                            ids = ' '.join( map(str, id) )
                        x = ll[0] + j * pitch[0]
                        y = ll[1] + i * pitch[1]
                        cmds.append( f"move volume {ids} midpoint location {x} {y} 0 except z" )
                    j = j + 1
                i = i + 1
        #FIXME rotate and tranlate
        return results

    #print( geom.root_universe )
    for cell in geom.root_universe._cells.values():
        process_node_or_fill( cell, w )

    if filename:
        with f as open( filename, "w" ):
            for x in cmds:
                f.write( x + "\n" )
            #f.write( "save sub5 \"geometry.cub\" overwrite" )

    if to_cubit:
        if cubit:
            cubit.cmd( "reset" )
            for x in cmds:
                cubit.cmd( x )

    #for x in cmds:
    #    print( x )

