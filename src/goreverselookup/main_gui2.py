# The main GUI class of GOReverseLookup.

import moderngl_window as mglw

class App(mglw.WindowConfig):
    window_size = (800, 600)
    resource_dir = '.'

    def key_event(self, key, action, modifiers):
        # Close the window when the ESC key is pressed
        if key == self.wnd.keys.ESCAPE:
            self.wnd.close()

    def render(self, time, frame_time):
        # Clear the screen with a solid color
        self.ctx.clear(0.2, 0.4, 0.7)

if __name__ == '__main__':
    mglw.run_window_config(App)