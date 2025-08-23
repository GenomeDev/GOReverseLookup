# The main GUI class of GOReverseLookup.

import moderngl_window as mglw
import numpy as np
import textwrap

vertex_shader=textwrap.dedent('''#version 330
in vec2 in_position;
void main() {
    gl_Position = vec4(in_position, 0.0, 1.0);
}
''')
        
fragment_shader=textwrap.dedent('''#version 330
out vec4 fragColor;
void main() {
    fragColor = vec4(0.0, 0.5, 1.0, 1.0);  # Some color
}
''')

vertices = np.array([
        -0.5, -0.5,
         0.5, -0.5,
         0.0,  0.5,
    ], dtype='f4')

class Window(mglw.WindowConfig):
    gl_version = (3, 3)
    window_size = (1920, 1080)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # Do initialization here
        self.prog = self.ctx.program(
            vertex_shader=vertex_shader,
            fragment_shader=fragment_shader
        )
        
        self.vao = self.ctx.vertex_array(
            self.prog,
            [(self.ctx.buffer(vertices).tobytes(), '2f', 'in_position')]
        )
        self.texture = self.ctx.texture(self.wnd.size, 4)

    def render(self, time, frametime):
        # This method is called every frame
        self.vao.render()
    
    def resize(self, width: int, height: int):
        print("Window was resized. buffer size is {} x {}".format(width, height))

    def close(self):
        print("The window is closing")

    def iconify(self, iconify: bool):
        print("Window was iconified:", iconify)
    
    def key_event(self, key, action, modifiers):
        # Key presses
        if action == self.wnd.keys.ACTION_PRESS:
            if key == self.wnd.keys.SPACE:
                print("SPACE key was pressed")

            # Using modifiers (shift and ctrl)

            if key == self.wnd.keys.Z and modifiers.shift:
                print("Shift + Z was pressed")

            if key == self.wnd.keys.Z and modifiers.ctrl:
                print("ctrl + Z was pressed")
            
            if key == self.wnd.keys.ESCAPE:
                print("ESCAPE key was pressed. Exiting.")
                self.wnd.close()
                       
        # Key releases
        elif action == self.wnd.keys.ACTION_RELEASE:
            if key == self.wnd.keys.SPACE:
                print("SPACE key was released")
    
    def mouse_position_event(self, x, y, dx, dy):
        print("Mouse position:", x, y, dx, dy)

    def mouse_drag_event(self, x, y, dx, dy):
        print("Mouse drag:", x, y, dx, dy)

    def mouse_scroll_event(self, x_offset: float, y_offset: float):
        print("Mouse wheel:", x_offset, y_offset)

    def mouse_press_event(self, x, y, button):
        print("Mouse button {} pressed at {}, {}".format(button, x, y))

    def mouse_release_event(self, x: int, y: int, button: int):
        print("Mouse button {} released at {}, {}".format(button, x, y))

    def unicode_char_entered(self, char: str):
        print('character entered:', char)

# Blocking call entering rendering/event loop
mglw.run_window_config(Window)


    