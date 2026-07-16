# Copyright 2026 Capytaine developers
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import logging
from rich.logging import RichHandler


def set_logging(level="INFO", force=False):
    """Configure logging with a nice Rich handler.

    If the root logger already has handlers (i.e., the user has set up their own
    logging configuration), this function does nothing unless ``force=True``.
    """
    if force or not logging.root.handlers:
        logging.basicConfig(
            level=level,
            format="%(message)s",
            handlers=[RichHandler(level=level, log_time_format="[%X]", show_path=False)],
            force=force,
        )
