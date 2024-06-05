import math
import sys
from openmc.region import Region, Complement, Intersection, Union
from openmc.surface import Halfspace
from openmc.lattice import Lattice, HexLattice

def flatten(S):
    if S == []:
        return S
    if isinstance(S[0], list):
        return flatten(S[0]) + flatten(S[1:])
    return S[:1] + flatten(S[1:])

def vector_to_euler_xyz(v):
    x, y, z = v
    phi = math.atan2(z, x)
    theta = math.acos(x / math.sqrt(x**2 + y**2))
    psi = math.atan2(y * math.cos(theta), x)

    # Ensure angles are within [0, 2*pi] range
    phi %= (2 * math.pi)
    theta %= (2 * math.pi)
    psi %= (2 * math.pi)

    oe = 180 / math.pi 
    return phi * oe, theta * oe, psi * oe

def to_cubit_journal(geom, seen=set(), world=[60,60,60], filename=None, to_cubit=False ):
    w = world
    cid = 1
    def lastid():
        nonlocal cid 
        id = cid
        cid = cid + 1
        return id
    cmds = []
    cmds.extend( [
        #"set echo off",
        "set graphics off",
        "set journal off",
        #"set undo off",
        #"set default autosizing off",   # especially when the CAD is complex (contains many spline surfaces) this can have a massive effect
        #"set info off",
        #"set warning off",
        ])
    def python_cmd( s ):
        cmds.append( s )
    def cubit_cmd( s ):
        cmds.append( s )
        #cmds.append( f'cubit.cmd( "{s}" )' )

    def emit_get_last_id( ):
        idn = lastid()
        ids = f"id{idn}"
        python_cmd( f'#{{ {ids} = Id("volume") }}' )
        return ids

    def rotate( id, x, y, z ):
        if nonzero( x, y, z ):
            phi, theta, psi = vector_to_euler_xyz( ( x, y, z ) )
            cubit_cmd( f"body {{ {id} }} rotate {phi} about Z" )
            cubit_cmd( f"body {{ {id} }} rotate {theta} about Y" )
            cubit_cmd( f"body {{ {id} }} rotate {psi} about X" )

    def nonzero(*args):
        return any(arg!= 0 for arg in args)

    def move( id, x, y, z ):
        if nonzero( x, y, z ):
           cubit_cmd( f"body {{ {id} }} move {x} {y} {z}" )

    def make_world_brick():
        pass
        
    def surface_to_cubit_journal(node, w, indent = 0, inner_world = None, hex = False ):
        def ind():
            return ' ' * (2*indent)
        if isinstance(node, Halfspace):
            #if not ( node.surface in seen ):
                seen.add( node.surface )
                surface = node.surface
               #print( ind(), "Surface:", surface )

                def reverse():
                    if node.side == '-':
                        #print( "REV", node, node.side )
                        return "reverse"
                    else:
                        #print( "NORMAL", node, node.side )
                        return ""

                if surface._type == "plane":
                    cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                    ids = emit_get_last_id()
                    phi, theta, psi = vector_to_euler_xyz( ( surface.coefficients['a'], surface.coefficients['b'], surface.coefficients['c'] ) )
                    cmds.append( f"body {{ { ids } }} rotate {phi} about Z" )
                    cmds.append( f"body {{ { ids } }} rotate {theta} about Y" )
                    cmds.append( f"body {{ { ids } }} rotate {psi} about X" )
                    cmds.append( f"body {{ { ids } }} move direction {surface.coefficients['a'] } { surface.coefficients['b']} {surface.coefficients['c']} distance {surface.coefficients['d']}" )
                    return ids
                elif surface._type == "x-plane":
                    cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                    ids = emit_get_last_id()
                    cmds.append( f"section volume {{ {ids} }} with xplane offset {surface.coefficients['x0']} {reverse()}")
                    #ids = emit_get_last_id()
                    return ids
                elif surface._type == "y-plane":
                    cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                    ids = emit_get_last_id()
                    cmds.append( f"section volume {{ {ids} }} with yplane offset {surface.coefficients['y0']} {reverse()}")
                    #ids = emit_get_last_id()
                    return ids
                elif surface._type == "z-plane":
                    cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                    ids = emit_get_last_id()
                    cmds.append( f"section volume {{ {ids} }} with zplane offset {surface.coefficients['z0']} {reverse()}")
                    #ids = emit_get_last_id()
                    return ids
                elif surface._type == "cylinder":
                    h = inner_world[2] if inner_world else w[2] 
                    cmds.append( f"cylinder height {h} radius {surface.coefficients['r']}")
                    ids = emit_get_last_id()
                    if node.side != '-':
                        wid = 0
                        if inner_world:
                            if hex:
                                cmds.append( f"create prism height {inner_world[2]} sides 6 radius { ( inner_world[0] / 2 ) }" )
                                wid = emit_get_last_id()
                                cmds.append( f"rotate vol {{ {wid} }} about z angle 30" )
                            else:
                                cmds.append( f"brick x {inner_world[0]} y {inner_world[1]} z {inner_world[2]}" )
                                wid = emit_get_last_id()
                        else:
                            cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                            wid = emit_get_last_id()
                        cmds.append( f"subtract vol {{ { ids } }} from vol {{ { wid } }}" )
                        rotate( wid, surface.coefficients['dx'], surface.coefficients['dy'], surface.coefficients['dz'] )
                        move( wid, surface.coefficients['x0'], surface.coefficients['y0'], surface.coefficients['z0'] )
                        return wid
                    rotate( ids, surface.coefficients['dx'], surface.coefficients['dy'], surface.coefficients['dz'] )
                    move( ids, surface.coefficients['x0'], surface.coefficients['y0'], surface.coefficients['z0'] )
                    return ids
                elif surface._type == "x-cylinder":
                    h = inner_world[0] if inner_world else w[0] 
                    cmds.append( f"cylinder height {h} radius {surface.coefficients['r']}")
                    ids = emit_get_last_id()
                    cmds.append( f"rotate volume {{ { ids } }} about y angle 90")
                    if node.side != '-':
                        wid = 0
                        if inner_world:
                            if hex:
                                cmds.append( f"create prism height {inner_world[2]} sides 6 radius { inner_world[0] / 2 } " )
                                wid = emit_get_last_id()
                                cmds.append( f"rotate vol {{ {wid} }} about z angle 30" )
                                cmds.append( f"rotate vol {{ {wid} }} about y angle 90")
                            else:
                                cmds.append( f"brick x {inner_world[0]} y {inner_world[1]} z {inner_world[2]}" )
                                wid = emit_get_last_id()
                        else:
                            cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                            wid = emit_get_last_id()
                        cmds.append( f"subtract vol {{ { ids } }} from vol {{ { wid } }}" )
                        move( wid, 0, surface.coefficients['y0'], surface.coefficients['z0'] )
                        return wid
                    move( ids, 0, surface.coefficients['y0'], surface.coefficients['z0'] )
                    return ids
                elif surface._type == "y-cylinder":
                    h = inner_world[1] if inner_world else w[1] 
                    cmds.append( f"cylinder height {h} radius {surface.coefficients['r']}")
                    ids = emit_get_last_id()
                    cmds.append( f"rotate volume {{ {ids} }} about x angle 90")
                    if node.side != '-':
                        wid = 0
                        if inner_world:
                            if hex:
                                cmds.append( f"create prism height {inner_world[2]} sides 6 radius { ( inner_world[0] / 2) }" )
                                wid = emit_get_last_id()
                                cmds.append( f"rotate vol {{ {wid} }} about z angle 30" )
                                cmds.append( f"rotate vol {{ {wid} }} about x angle 90")
                            else:
                                cmds.append( f"brick x {inner_world[0]} y {inner_world[1]} z {inner_world[2]}" )
                                wid = emit_get_last_id()
                        else:
                            cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                            wid = emit_get_last_id()
                        cmds.append( f"subtract vol {{ {id} }} from vol {{ {wid} }}" )
                        move( wid, surface.coefficients['x0'], 0, surface.coefficients['z0'] )
                        return wid
                    move( ids, surface.coefficients['x0'], 0, surface.coefficients['z0'] )
                    return ids
                elif surface._type == "z-cylinder":
                    h = inner_world[2] if inner_world else w[2] 
                    cmds.append( f"cylinder height {h} radius {surface.coefficients['r']}")
                    ids = emit_get_last_id()
                    if node.side != '-':
                        wid = 0
                        if inner_world:
                            if hex:
                                cmds.append( f"create prism height {inner_world[2]} sides 6 radius { ( inner_world[0] / 2 ) }" )
                                wid = emit_get_last_id()
                                cmds.append( f"rotate vol {{ {wid} }} about z angle 30" )
                            else:
                                cmds.append( f"brick x {inner_world[0]} y {inner_world[1]} z {inner_world[2]}" )
                                wid = emit_get_last_id()
                        else:
                            cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                            wid = emit_get_last_id()
                        cmds.append( f"subtract vol {{ { ids } }} from vol {{ { wid } }}" )
                        move( wid, surface.coefficients['x0'], surface.coefficients['y0'], 0 )
                        return wid
                    move( ids, surface.coefficients['x0'], surface.coefficients['y0'], 0 )
                    return ids
                elif surface._type == "sphere":
                    cmds.append( f"sphere redius {surface.coefficients['r']}")
                    ids = emit_get_last_id()
                    move( ids, surface.coefficients['x0'], surface.coefficients['y0'], surface.coefficients['zy0'] )
                    pass
                elif surface._type == "cone":
                    raise "cone not implemented"
                    pass
                elif surface._type == "x-cone":
                    cmds.append( f"create frustum height {w[0]} radius {math.sqrt(surface.coefficients['r2']*w[0])} top 0")
                    ids = emit_get_last_id()
                    cmds.append( f"rotate volume {{ {ids} }} about y angle 90")
                    if node.side != '-':
                        cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                        wid = emit_get_last_id()
                        cmds.append( f"subtract vol {{ {ids} }} from vol {{ {wid} }}" )
                        move( wid, surface.coefficients['x0'], surface.coefficients['y0'], surface.coefficients['z0'] )
                        return wid
                    move( ids, surface.coefficients['x0'], surface.coefficients['y0'], surface.coefficients['z0'] )
                    return ids
                elif surface._type == "y-cone":
                    cmds.append( f"create frustum height {w[1]} radius {math.sqrt(surface.coefficients['r2']*w[1])} top 0")
                    ids = emit_get_last_id()
                    cmds.append( f"rotate volume {{ {ids} }} about x angle 90")
                    if node.side != '-':
                        cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                        wid = emit_get_last_id()
                        cmds.append( f"subtract vol {{ {ids} }} from vol {{ {wid} }}" )
                        move( wid, surface.coefficients['x0'], surface.coefficients['y0'], surface.coefficients['z0'] )
                        return wid
                    move( ids, surface.coefficients['x0'], surface.coefficients['y0'], surface.coefficients['z0'] )
                    return ids
                elif surface._type == "z-cone":
                    cmds.append( f"create frustum height {w[2]} radius {math.sqrt(surface.coefficients['r2']*w[2])} top 0")
                    ids = emit_get_last_id()
                    if node.side != '-':
                        cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                        wid = emit_get_last_id()
                        cmds.append( f"subtract vol {{ {ids} }} from vol {{ {wid} }}" )
                        move( wid, surface.coefficients['x0'], surface.coefficients['y0'], surface.coefficients['z0'] )
                        return wid
                    move( ids, surface.coefficients['x0'], surface.coefficients['y0'], surface.coefficients['z0'] )
                    return ids
                elif surface._type == "x-torus":
                    cmds.append( f"torus major radius {surface.coefficients['r']} minor radius {surface.coefficients['r']}")
                    ids = emit_get_last_id()
                    cmds.append( f"rotate volume {{ {ids} }} about y angle 90")
                    if node.side != '-':
                        cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                        wid = emit_get_last_id()
                        cmds.append( f"subtract vol {{ {ids} }} from vol {{ {wid} }}" )
                        move( wid, surface.coefficients['x0'], surface.coefficients['y0'], surface.coefficients['z0'] )
                        return wid
                    move( ids, surface.coefficients['x0'], surface.coefficients['y0'], surface.coefficients['z0'] )
                    return ids
                elif surface._type == "y-torus":
                    cmds.append( f"torus major radius {surface.coefficients['r']} minor radius {surface.coefficients['r']}")
                    ids = emit_get_last_id()
                    cmds.append( f"rotate volume {{ {ids} }} about x angle 90")
                    if node.side != '-':
                        cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                        wid = emit_get_last_id()
                        cmds.append( f"subtract vol {{ {id} }} from vol {{ {wid} }}" )
                        move( wid, surface.coefficients['x0'], surface.coefficients['y0'], surface.coefficients['z0'] )
                        return wid
                    return ids
                elif surface._type == "z-torus":
                    cmds.append( f"torus major radius {surface.coefficients['r']} minor radius {surface.coefficients['r']}")
                    ids = emit_get_last_id()
                    if node.side != '-':
                        cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
                        wid = emit_get_last_id()
                        cmds.append( f"subtract vol {{ {ids} }} from vol {{ {wid} }}" )
                        move( wid, surface.coefficients['x0'], surface.coefficients['y0'], surface.coefficients['z0'] )
                        return wid
                    move( ids, surface.coefficients['x0'], surface.coefficients['y0'], surface.coefficients['z0'] )
                    return ids
                else:
                    raise f"{surface.type} not implemented"
            #else:
            #    print( ind(), node )

        elif isinstance(node, Complement):
            print( "Complement:" )
            id = surface_to_cubit_journal(node.node, w, indent + 1, inner_world )
            cmds.append( f"brick x {w[0]} y {w[1]} z {w[2]}" )
            wid = emit_get_last_id()
            cmds.append( f"subtract vol {{ {id} }} from vol {{ {wid} }}" )
            return emit_get_last_id()
        elif isinstance(node, Intersection):
            #print( ind(), "Intersection:" )
            surfaces = []
            for subnode in node:
                #print( ind(), "Subnode:", subnode )
                id = surface_to_cubit_journal( subnode, w, indent + 1, inner_world,)
                surfaces.append( id )
            if len( surfaces) > 1:
                last = surfaces[0]
                for s in surfaces[1:]:
                    cmds.append( f"intersect {{ {last} }} {{ {s} }}" )
                    last = s
                if inner_world:
                    cmds.append( f"brick x {inner_world[0]} y {inner_world[1]} z {inner_world[2]}" )
                    iwid = emit_get_last_id()
                    cmds.append( f"intersect {{ {last} }} {{ {iwid} }}" )
                    return emit_get_last_id()
            return surfaces[-1]
        elif isinstance(node, Union):
            #print( ind(), "Union:" )
            surfaces = []
            for subnode in node:
                ##print( ind(), "Subnode:", subnode )
                id = surface_to_cubit_journal( subnode, w, indent + 1, inner_world )
                surfaces.append( id )
            if len( surfaces) > 1:
                last = surfaces[0]
                for s in surfaces[1:]:
                    cmds.append( f"unite {{ {last} }} {{ {s} }}" )
                    last = s
                if inner_world:
                    cmds.append( f"brick x {inner_world[0]} y {inner_world[1]} z {inner_world[2]}" )
                    iwid = emit_get_last_id()
                    cmds.append( f"intersect {{ {last} }} {{ {iwid} }}" )
                    return emit_get_last_id()
            return surfaces[-1]
        elif isinstance(node, None):
            pass
        else:
            #print( node )
            ##print( ind(), "Composite:", node )
            for subnode in node:
                ##print( ind(), "Subnode:", subnode )
                surface_to_cubit_journal( subnode, w, indent + 1, inner_world )

    def process_node_or_fill( node, w, indent = 0, offset = [0, 0], inner_world = None, outer_ll = None, hex = False ):
        def ind():
            return ' ' * (2*indent)
        #if not ( node in seen ):
        #print(  ind(), "Node:", node )
        seen.add( node )
        results = []
        if hasattr( node, "region" ) and not ( hasattr( node, "fill" ) and isinstance(node.fill, Lattice) ):
            if node.region != None:
                #print(  ind(), "Region:", node.region )
                id = surface_to_cubit_journal( node.region, w, indent, inner_world, hex = hex )
                results.append( id )
            elif hex:
                cmds.append( f"create prism height {inner_world[2]} sides 6 radius { ( inner_world[0] / 2) }" )
                wid = emit_get_last_id()
                cmds.append( f"rotate vol {{ {wid} }} about z angle 30" )
                results.append( wid )

        if hasattr( node, "fill" ) and isinstance(node.fill, Lattice):
            #print(  ind(), "Fill:", node.fill )
            id = process_node_or_fill( node.fill, w, indent + 1, offset, inner_world )
            results.append( id )
        if hasattr( node, "universes" ):
            #print(  ind(), "Universes:", node )
            pitch = node._pitch

            if isinstance( node, HexLattice ) :
                three_d_hex_lattice = len( node._pitch ) > 1
                #three_d_hex_lattice = len( node.center ) > 2
                #3d hex lattice
                if three_d_hex_lattice:
                    center = [ node.center[0], node.center[1], node.center[1] ]
                    # ll = [ node.lower_left[0], node.lower_left[1] ]
                    # if outer_ll:
                    #     ll = outer_ll
                    # ll[0] = ll[0] + pitch[0] / 2
                    # ll[1] = ll[1] + pitch[1] / 2
                    ii = 0
                    for uss in node.universes:
                        z = ii * pitch[1]
                        j = 0
                        for us in uss:
                            k = 0
                            ring_id = ( len( uss ) - j -1 )
                            #print( "RING_ID", ring_id )
                            def draw_hex_cell( n, cell, x, y ):
                                #print( i, j, k, len( node.universes ), len( uss), len( us ), x, y )
                                id = process_node_or_fill( cell, [ w[0], w[1], w[2] ], indent + 1, offset = [x,y,z], inner_world=[ pitch[0], pitch[0], pitch[1] ], outer_ll=outer_ll if outer_ll else [ node.center[0], node.center[1] ], hex = True )
                                ids = str( id )
                                if isinstance( id, list ):
                                    ids = ' '.join( map(str, id) )
                                else:
                                    pass
                                if ids != '':
                                    cmds.append( f"move volume {{ {ids} }} midpoint location {x} {y} {z}" )
                            side_to_side_diameter =  pitch[0]/2 * math.sqrt( 3 )
                            center_to_mid_side_diameter = ( ( pitch[0] / 2 ) * math.sin( math.pi / 6 ) ) + pitch[0] / 2
                            #print( "diameter", pitch[0] )
                            #print( "side_to_side_diameter", side_to_side_diameter )
                            #print( "center_to_mid_side_diameter", center_to_mid_side_diameter)
                            if ring_id < 2:
                                for u in us:
                                    for n, cell in u._cells.items():
                                        #print( n, cell )
                                        theta = 2 * math.pi * -k / len( us ) + math.pi / 2
                                        r = ( len( uss ) - j -1 ) * side_to_side_diameter
                                        x = r * math.cos( theta )
                                        y = r * math.sin( theta )
                                        draw_hex_cell( n, cell, x, y )
                                    k = k + 1
                            else:
                                x = 0
                                x_pos = 0
                                y_pos = 0
                                r = ring_id
                                for i in range( r, 0, -1 ):
                                    x_pos = x * center_to_mid_side_diameter; 
                                    y_pos = ring_id * side_to_side_diameter - ( x ) * 0.5 * side_to_side_diameter; 
                                    for n, cell in us[k]._cells.items():
                                        draw_hex_cell( n, cell, x_pos, y_pos )
                                    #print( r, k, x, x_pos, y_pos )
                                    k = k + 1
                                    x = x + 1
                                y_pos = ring_id * side_to_side_diameter - ( x ) * 0.5 * side_to_side_diameter; 
                                for i in range( r, 0, -1 ):
                                    x_pos = x * center_to_mid_side_diameter; 
                                    for n, cell in us[k]._cells.items():
                                        draw_hex_cell( n, cell, x_pos, y_pos )
                                    #print( r, k, x, x_pos, y_pos )
                                    y_pos = y_pos - side_to_side_diameter; 
                                    k = k + 1
                                for i in range( r, 0, -1 ):
                                    x_pos = x * center_to_mid_side_diameter; 
                                    y_pos = - ring_id * side_to_side_diameter + ( x ) * 0.5 * side_to_side_diameter; 
                                    for n, cell in us[k]._cells.items():
                                        draw_hex_cell( n, cell, x_pos, y_pos )
                                    #print( r, k, x, x_pos, y_pos )
                                    k = k + 1
                                    x = x - 1
                                for i in range( r, 0, -1 ):
                                    x_pos = x * center_to_mid_side_diameter; 
                                    y_pos = - ring_id * side_to_side_diameter - ( x ) * 0.5 * side_to_side_diameter; 
                                    for n, cell in us[k]._cells.items():
                                        draw_hex_cell( n, cell, x_pos, y_pos )
                                    #print( r, k, x, x_pos, y_pos )
                                    k = k + 1
                                    x = x - 1
                                y_pos = - ring_id * side_to_side_diameter - ( x ) * 0.5 * side_to_side_diameter; 
                                for i in range( r, 0, -1 ):
                                    x_pos = x * center_to_mid_side_diameter; 
                                    for n, cell in us[k]._cells.items():
                                        draw_hex_cell( n, cell, x_pos, y_pos )
                                    #print( r, k, x, x_pos, y_pos )
                                    y_pos = y_pos + side_to_side_diameter; 
                                    k = k + 1
                                for i in range( r, 0, -1 ):
                                    x_pos = x * center_to_mid_side_diameter; 
                                    y_pos = ring_id * side_to_side_diameter + ( x ) * 0.5 * side_to_side_diameter; 
                                    for n, cell in us[k]._cells.items():
                                        draw_hex_cell( n, cell, x_pos, y_pos )
                                    #print( r, k, x, x_pos, y_pos )
                                    k = k + 1
                                    x = x + 1
                            j = j + 1
                        ii = ii + 1
                #2d hex lattice
                else:
                    #print( node )
                    center = [ node.center[0], node.center[1] ]
                    i = 0
                    for us in node.universes:
                        j = 0
                        for u in us:
                            for n, cell in u._cells.items():
                                    #print( n, cell )
                                    theta = 2 * math.pi * -j / len( us ) + math.pi / 2
                                    r = ( len( uss ) - i -1 ) * pitch[0]
                                    x = r * math.cos( theta )
                                    y = r * math.sin( theta )
                                    #print( n, i, j, k, len( node.universes ), len( uss), len( us ), x, y, theta )
                                    id = process_node_or_fill( cell, [ w[0], w[1], w[2] ], indent + 1, offset = [x,y], inner_world=[ pitch[0], pitch[0], pitch[1] ], outer_ll=outer_ll if outer_ll else [ node.center[0], node.center[1] ], hex = True )
                                    ids = str( id )
                                    if isinstance( id, list ):
                                        ids = ' '.join( map(str, id) )
                                        #results.extend( id )
                                    else:
                                        #results.append( id )
                                        pass
                                    if ids != '':
                                        cmds.append( f"move volume {{ {ids} }} midpoint location {x} {y} {z}" )
                            j = j + 1
                        i = i + 1

            elif isinstance( node, Lattice ):
                ll = [ node.lower_left[0], node.lower_left[1] ]
                if outer_ll:
                    ll = outer_ll
                ll[0] = ll[0] + pitch[0] / 2
                ll[1] = ll[1] + pitch[1] / 2
                i = 0
                for us in node.universes:
                    j = 0
                    for u in us:
                        for n, cell in u._cells.items():
                            x = ll[0] + j * pitch[0] + offset[0]
                            y = ll[1] + i * pitch[1] + offset[1]
                            #print(  ind(), "UCell:", n, cell )
                            id = process_node_or_fill( cell, [ w[0], w[1], w[2] ], indent + 1, offset = [x,y], inner_world=[ pitch[0], pitch[1], w[2] ], outer_ll=outer_ll if outer_ll else [ node.lower_left[0], node.lower_left[1] ] )
                            ids = str( id )
                            if isinstance( id, list ):
                                ids = ' '.join( map(str, id) )
                                #results.extend( id )
                            else:
                                #results.append( id )
                                pass
                            if ids != '':
                                cmds.append( f"move volume {{ {ids} }} midpoint location {x} {y} 0 except z" )
                        j = j + 1
                    i = i + 1
        #FIXME rotate and tranlate
        r = flatten( results )
        if len( r ) > 0  and node.name:
             cmds.append( f"volume {{ {r[0]} }} name \"{node.name}\"" )
        #print( r )
        return r

    #print( geom.root_universe )
    for cell in geom.root_universe._cells.values():
        process_node_or_fill( cell, w )

    if filename:
        with open( filename, "w" ) as f:
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

